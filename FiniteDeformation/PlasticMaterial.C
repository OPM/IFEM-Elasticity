// $Id$
//==============================================================================
//!
//! \file PlasticMaterial.C
//!
//! \date Mar 16 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Elasto-plastic material models.
//!
//==============================================================================

#include "PlasticMaterial.h"
#include "FiniteElement.h"
#include "TimeDomain.h"
#include "Function.h"
#include "Vec3Oper.h"
#include "IFEM.h"


PlasticMaterial::PlasticMaterial (const RealArray& p, const ScalarFunc* hcurve)
  : pMAT(p), hardening(hcurve), iAmIntegrating(false), firstStep(true), iP2(0)
{
  if (pMAT.size() < 11) pMAT.resize(11,0.0);

  double Emod = pMAT[0];
  double nu   = pMAT[1];
  double Bmod, Smod;
  if (nu > 0.5)
  {
    // Calculate E and nu from the Bulk and Shear moduli
    Bmod = Emod;
    Smod = nu;
    Emod = 9.0*Bmod*Smod/(3.0*Bmod + Smod);
    nu   = (1.5*Bmod - Smod)/(3.0*Bmod + Smod);
  }
  else if (nu < 0.5)
  {
    // Calculate the Bulk and Shear moduli from E and nu
    Bmod = Emod / (3.0 - 6.0*nu);
    Smod = Emod / (2.0 + 2.0*nu);
  }
  else
  {
    Bmod = 0.0;
    Smod = Emod / 3.0;
  }

  pMAT[0] = Emod;
  pMAT[1] = nu;
  pMAT.insert(pMAT.begin()+4,Bmod);
  pMAT.insert(pMAT.begin()+5,Smod);
}


PlasticMaterial::~PlasticMaterial ()
{
  for (ResultPoint* pt : itgPoints) delete pt;
  for (ResultPoint* pt : resPoints) delete pt;
  delete hardening;
}


void PlasticMaterial::printLog () const
{
  IFEM::cout <<"PlasticMaterial: pMAT =";
  for (double v : pMAT) IFEM::cout <<" "<< v;
  IFEM::cout << std::endl;
}


void PlasticMaterial::initIntegration (size_t nGP)
{
  itgPoints.resize(nGP,nullptr);
}


void PlasticMaterial::initIntegration (const TimeDomain& prm)
{
#if INT_DEBUG > 0
  std::cout <<"PlasticMaterial::initIntegration: "<< itgPoints.size()
            << std::endl;
#endif

  iAmIntegrating = true;
  firstStep = prm.first;

#if INT_DEBUG > 0
  int nUpdated = 0;
  for (ResultPoint* itgPt : itgPoints)
    if (itgPt && itgPt->updateState(!prm.first && prm.it == 0))
      nUpdated++;
  std::cout <<"PlasticMaterial::initIntegration: History updated "<< nUpdated
            << std::endl;
#else
  for (ResultPoint* itgPt : itgPoints)
    if (itgPt) itgPt->updateState(!prm.first && prm.it == 0);
#endif
}


void PlasticMaterial::initResultPoints ()
{
#if INT_DEBUG > 0
  std::cout <<"PlasticMaterial::initResultPoints: "<< iP2 << std::endl;
#endif

  iP2 = 0;
  iAmIntegrating = false;

#if INT_DEBUG > 0
  int nUpdated = 0;
  for (ResultPoint* resPt : resPoints)
    if (resPt->updateState())
      nUpdated++;
  std::cout <<"PlasticMaterial::initResultPoints: History updated "<< nUpdated
            << std::endl;
#else
  for (ResultPoint* resPt : resPoints)
    resPt->updateState();
#endif
}


bool PlasticMaterial::evaluate (Matrix& C, SymmTensor& sigma, double& U,
                                const FiniteElement& fe, const Vec3&,
                                const Tensor& F, const SymmTensor& eps,char iop,
                                const TimeDomain* prm, const Tensor* Fpf) const
{
  C.resize(sigma.size(),sigma.size());

  bool ok = true;
  if (iAmIntegrating)
  {
    if (!prm)
    {
      static TimeDomain dummy(1);
      dummy.first = firstStep;
      prm = &dummy;
    }
    size_t iP1 = fe.iGP;
    if (iP1 >= itgPoints.size())
    {
      std::cerr <<" *** PlasticMaterial::evaluate: Integration point "<< iP1+1
                <<" out of range [1,"<< itgPoints.size() <<"]."<< std::endl;
      return false;
    }
    else if (!itgPoints[iP1])
    {
      if (prm->it == 0 && prm->first)
        itgPoints[iP1] = new ResultPoint(this,F.dim());
      else
      {
        std::cerr <<" *** PlasticMaterial::evaluate: Integration point "<< iP1
                  <<" does not exist."<< std::endl;
        return false;
      }
    }

    if (prm->it == 0 && !prm->first)
      itgPoints[iP1]->Fp = F;

#if INT_DEBUG > 0
    std::cout <<"PlasticMaterial: Evaluating itg.point #"<< iP1+1 << std::endl;
#endif
    if (itgPoints[iP1]->evaluate(C,sigma,F,*prm) != 0)
      ok = false;

    // Calculate principal stresses, etc.
    itgPoints[iP1]->principalStress(sigma);
  }
  else // Result evaluation
  {
    while (resPoints.size() <= iP2)
      resPoints.push_back(new ResultPoint(this,F.dim()));

#if INT_DEBUG > 0
    std::cout <<"PlasticMaterial: Evaluating result point #"<< iP2 << std::endl;
#endif
    // Always invoke with iter = 1 in result evaluation
    if (resPoints[iP2]->evaluate(C,sigma,F,TimeDomain(1,false)) < 0)
      ok = false; // no error return here on material point divergence

    // Calculate principal stresses, etc.
    resPoints[iP2]->principalStress(sigma);

    // Assume only one evaluation per increment; always update Fp
    resPoints[iP2++]->Fp = F;
  }

  if (iop > 1)
  {
    // Transform to 2nd Piola-Kirchhoff stresses,
    // via pull-back to reference configuration
    Tensor Fi(Fpf ? *Fpf : F);
    double J = Fi.inverse();
    if (iop == 2)
    {
      sigma.transform(Fi); // sigma = F^-1 * sigma * F^-t
      sigma *= J;
      //TODO: Also pull-back the C-matrix (Total Lagrange formulation)
      std::cerr <<" *** PlasticMaterial::evaluate: Not available for"
                <<" Total Lagrangian formulation, sorry."<< std::endl;
      ok = false;
    }
    else
    {
      SymmTensor S(sigma); // sigma should be Cauchy stress when iop=3
      S.transform(Fi);     // S = F^-1 * sigma * F^-t
      if (iAmIntegrating)
        U = itgPoints[fe.iGP]->energyIntegral(S,eps)*J;
      else
        U = resPoints[iP2-1]->energyIntegral(S,eps)*J;
    }
  }

  return ok;
}


bool PlasticMaterial::diverged (size_t iP1) const
{
  if (iP1 > 0 && --iP1 < itgPoints.size())
    return itgPoints[iP1] ? itgPoints[iP1]->diverged() : false;
  else for (ResultPoint* itgPt : itgPoints)
    if (itgPt && itgPt->diverged())
      return true;

  return false;
}


double PlasticMaterial::getInternalVar (int idx, char* label, size_t iP1) const
{
  ResultPoint* p = nullptr;
  if (iAmIntegrating)
  {
    if (iP1 < itgPoints.size())
      p = itgPoints[iP1];
  }
  else if (iP2 > 0)
    if (iP2 <= resPoints.size())
      p = resPoints[iP2-1];

  switch (idx) {
  case 1:
    if (label) strcpy(label,"E_pp"); // Equivalent plastic strain
    return p ? (p->diverged() ? 0.0 : p->getVariable(6)) : 0.0;
  case 2:
    if (label) strcpy(label,"s_h"); // Mean stress
    return p ? p->getMeanStress() : 0.0;
  case 3:
  case 4:
  case 5:
    if (label) sprintf(label,"s_%d",idx-2); // Principal stress
    return p ? p->getPrnStress(idx-3) : 0.0;
  case 6:
    if (label) strcpy(label,"T"); // Stress triaxiality
    return p ? p->getTriax() : 0.0;
  case 7:
    if (label) strcpy(label,"L"); // Lode parameter
    return p ? p->getLode() : 0.0;
  default:
    if (label) strcpy(label,"zero");
    return 0.0;
  }
}


PlasticMaterial::PlasticPoint::PlasticPoint (const PlasticMaterial* prm,
                                             unsigned short int n)
  : pMAT(prm->pMAT), hfn(prm->hardening), updated(false), Ep(0), Sp(0), Fp(n)
{
  // Initialize the history variables
  HVc.fill(0.0);
  HVp.fill(0.0);
  Fp = HVp[0] = HVp[1] = HVp[2] = 1.0;
  Up = 0.0;
}


bool PlasticMaterial::PlasticPoint::updateState (bool updateVars)
{
  if (!updated) return false;
  updated = false;

  // Update history variables with values of the new converged solution
  if (updateVars)
    HVp = HVc;
  else
    return false;

#if INT_DEBUG > 0
  std::cout <<"PlasticMaterial: converged HV =";
  for (double v : HVp) std::cout <<" "<< v;
  std::cout << std::endl;
#endif

  return true;
}


int PlasticMaterial::PlasticPoint::evaluate (Matrix& C, SymmTensor& sigma,
                                             const Tensor& Fc,
                                             const TimeDomain& prm) const
{
  const int   itmax = 50;
  const double tolb = 1.0e-8;
  const double tolc = 1.0e-9;

  const double one2  = 0.5;
  const double one3  = 1.0/3.0;
  const double two3  = 2.0/3.0;
  const double sqt23 = sqrt(two3);

  const std::array<int,6> p1 = { 1,2,3,1,2,3 };
  const std::array<int,6> p2 = { 1,2,3,2,3,1 };

  HistoryVars& HVupd = const_cast<PlasticPoint*>(this)->HVc;
  double& Epp = HVupd[6]; // Accumulated plastic strain

  // Restore history variables from the previous, converged configuration
  HVupd = HVp;

#if INT_DEBUG > 0
  std::streamsize oldPrec = std::cout.precision(8);
  std::cout <<"PlasticMaterial: iter="<< prm.it <<" first="<< (int)prm.first;
  std::cout <<"\nPlasticMaterial::Fc =\n"<< Fc;
  std::cout <<"PlasticMaterial::Fp =\n"<< Fp;
  std::cout <<"PlasticMaterial::HV(in)  =";
  for (double v : HVc) std::cout <<" "<< v;
  std::cout << std::endl;
#endif

  // Fetch some material parameters

  const double Bmod = pMAT[4]; // Bulk modulus
  const double Smod = pMAT[5]; // Shear modulus
  const double Y0   = sqt23 * pMAT[9]; // Radius of yield = sqrt(2/3)*Sig_y
  const double Hk   = two3  * pMAT[7]; // Kinematic hardening = 2/3 * H_kin
  const double Hkr  = Hk > 0.0 ? 1.0 / Hk : 0.0;

  // Check state for iterations

  bool state = true;
  double dyld = 0.0;
  int istart = static_cast<int>(pMAT[12]); // Start state
  if (prm.it == 0)
  {
    // First iteration in step
    if (istart == 0) // Elastic state requested
      state = false;
    else             // Previous state requested
      dyld = 1.0e-8 * Y0;
  }
  if (prm.first && istart > 1 && prm.it <= istart)
    state = false;

  // Calculate the inverse of the previous deformation gradient
  // and the determinant of the current one

  Tensor Fi(Fp);
  double Jp = Fi.inverse();
  double Jc = Fc.det();
  if (Jp == 0.0 || Jc == 0.0)
  {
    std::cerr <<" *** PlasticMaterial::evaluate: "
              <<" Singular/zero deformation gradient(s)\n"<< Fp << Fc;
    return -999;
  }

  // Calculate the elastic left Cauchy-Green tensor in current configuration:
  // be = (Fc*Fi) * be * (Fc*Fi)^T

  SymmTensor be(RealArray(HVc.data(),HVc.data()+6));
  be.transform(Fi.preMult(Fc));
#if INT_DEBUG > 1
  std::cout <<"\nPlasticMaterial: be =\n"<< be;
#endif

  // Calculate principal stretches and directions

  Vec3 ll2_tr;
  Tensor nn_tr(3);
  if (!be.principal(ll2_tr,nn_tr))
    return -1;

  // Calculate trial Kirchhoff stress (pressure and deviatoric)

  Vec3 eps_tr;
  for (int i = 0; i < 3; i++)
    eps_tr[i] = log(sqrt(ll2_tr[i]));
  double th_tr = eps_tr.sum();
  double vol_tr = th_tr * one3;

#if INT_DEBUG > 1
  std::cout <<"\nPlasticMaterial: ll2_tr = "<< ll2_tr
            <<"\n                 eps_tr = "<< eps_tr <<" th_tr = "<< th_tr
            <<"\nPlasticMaterial: nn_tr\n"<< nn_tr;
#endif

  double Eppn  = Epp;
  double pp    = Bmod * th_tr; // Pressure: K*th_tr
  Vec3   tt    = 2.0*Smod * (eps_tr - vol_tr); // Trial deviatoric stress
  Vec3   tau   = tt + pp;
  Vec3   alp   = Hk * Vec3(HVc.data()+7);
  Vec3   alp_n = alp;
  SymmTensor dtde(3);

  // Deviatoric: ta = tt - alp_dev

  double aatr = alp.sum()*one3;

  Vec3 ta = tt - alp + aatr; // Trial Sigma = ss - alp

  // Compute stress invariant and yield function

  double f1, f2, f3, f11, f22, f33, f12, f13, f23, Ypr, YY, yield;
  double I1 = 3.0 * (pp-aatr);
  double J2 = ta.length2() * one2;
  double J3 = (ta.x*ta.x*ta.x + ta.y*ta.y*ta.y + ta.z*ta.z*ta.z) * one3;
  if (!this->yfunc(false,Epp,I1,J2,J3,
                   f1,f2,f3,f11,f22,f33,f12,f13,f23,Ypr,YY,yield))
    return -2;

  int i, j, a, b, ierr = 0;

  // Check yield

  if (yield > -dyld && state)
  {
    // Plastic step --> return map

    const double K3inv = one3 / Bmod;
    const double G2inv = one2 / Smod;
    const double d_el  = (K3inv - G2inv) * one3;

    double vol_e = K3inv * pp;
    Vec3 ee_e  = G2inv * tt;
    Vec3 eps_e = ee_e + vol_e;
    Vec3 nn, aa, bb, r1, r2;

    double epse = eps_tr.asum();
    double rerr = 1.0 + tolc;
    double gam  = 0.0;

    const size_t nvar = Hk > 0.0 ? 7 : 4;
    std::vector<int> iPivot(nvar);
    Matrix tres(nvar,nvar);
    Vector res(nvar);
    SymmTensor fss(3);
    Vec3 ll2;

    for (int it = 1; it < itmax && rerr >= tolc; it++)
    {
      if (!this->yfunc(true,Epp,I1,J2,J3,
                       f1,f2,f3,f11,f22,f33,f12,f13,f23,Ypr,YY,yield))
        return -2;

      double xx = f1 - f3*two3*J2;
      for (i = 0; i < 3; i++)
        nn[i] = xx + (f2 + f3*ta[i])*ta[i];

      r1 = eps_e - eps_tr + gam*nn;
      r2 = (alp_n - alp)*Hkr + gam*nn;
      rerr = r1.asum()/epse + (yield < 0.0 ? -yield : yield)/Y0;

      // Construct local tangent matrix

      double f112 = f11 - one3*f2;
      double f123 = f12 - two3*f3;

      for (i = 0; i < 3; i++)
      {
        aa[i] =   f23*ta[i] + f13;
        bb[i] = ta[i]*ta[i] - two3*J2;
      }

      for (a = 1; a <= 3; a++)
      {
        for (int b = a; b <= 3; b++)
          fss(a,b) = (f112 +
                      f22 * ta(a) * ta(b) +
                      f33 * bb(a) * bb(b) +
                      f123 * (ta(b) + ta(a)) +
                      aa(a) * bb(b) + bb(a) * aa(b)) * gam;
        fss(a,a) += (f2 + 2.0*f3*ta(a)) * gam;
      }

      for (a = 1; a <= 3; a++)
      {
        for (b = 1; b <= 3; b++)
          tres(a,b) = fss(a,b) + d_el;
        tres(a,a) += G2inv;
      }

      if (Hk > 0.0) // Kinematic and isotropic hardening
        for (a = 1; a <= 3; a++)
        {
          for (b = 1; b <= 3; b++)
          {
            tres(a  ,b+3) = tres(a+3,b) = -fss(a,b);
            tres(a+3,b+3) = fss(a,b);
          }

          tres(a+3,a+3) += 1.0/Hk;
          tres(a  ,7) = tres(7,a  ) =  nn(a);
          tres(a+3,7) = tres(7,a+3) = -nn(a);
          res(a)      = -r1(a);
          res(a+3)    = -r2(a);
        }

      else // Isotropic hardening only
        for (a = 1; a <= 3; a++)
        {
          tres(a,4) = tres(4,a) = nn(a);
          res(a)    = -r1(a);
        }

      tres(nvar,nvar) = -Ypr;
      res(nvar)       = -yield;
#if INT_DEBUG > 2
      std::cout <<"\nPlasticMaterial: tres"<< tres
                <<"\nPlasticMaterial: res"<< res;
#endif
      iPivot.front() = 0;
      if (!utl::solve(tres,res,&iPivot))
        return -3;

#if INT_DEBUG > 2
      std::cout <<"\nPlasticMaterial: tau = "<< tau
                <<"\n\nPlasticMaterial: dsol"<< res;
#endif
      tau += Vec3(res.ptr());
      gam += res[nvar-1];

      // Accumulated plastic strain
      Epp = Eppn + sqt23*gam;

#if INT_DEBUG > 2
      std::cout <<"\nPlasticMaterial: tau = "<< tau
                <<"  gam = "<< gam <<"  Epp = "<< Epp << std::endl;
#endif

      // Back stress
      if (nvar > 4)
        alp += Vec3(res.ptr()+3);

      // Volumetric-deviatoric Kirchhoff stress and stress invariants

      pp = tau.sum()*one3;
      tt = tau - pp;

      aatr = alp.sum()*one3;
      ta = tt - alp + aatr;

      I1 = 3.0 * (pp-aatr);
      J2 = ta.length2() * one2;
      J3 = (ta.x*ta.x*ta.x + ta.y*ta.y*ta.y + ta.z*ta.z*ta.z) * one3;

      // Volumetric-deviatoric logarithmic strain

      vol_e = K3inv * pp;
      ee_e  = G2inv * tt;
      eps_e = ee_e  + vol_e;
    }

    if (rerr >= tolc && prm.it > 0)
    {
      std::cout <<"  ** Warning: No convergence in plasticity iterations, "
                << rerr <<" " << tolc <<"\n     iter = "<< prm.it
                <<"\n     Epp  = "<< Epp <<"\n     Ypr  = "<< Ypr << std::endl;
      ierr = prm.it;
    }

    // Elastic left Cauchy-Green tensor and plastic accumulated strain

    for (a = 1; a <= 3; a++)
      ll2(a) = exp(2.0*eps_e(a));

    double* Ecg = HVupd.data();
    for (i = 0; i < 6; i++)
      Ecg[i] = (ll2.x * nn_tr(p1[i],1) * nn_tr(p2[i],1) +
                ll2.y * nn_tr(p1[i],2) * nn_tr(p2[i],2) +
                ll2.z * nn_tr(p1[i],3) * nn_tr(p2[i],3));

    // Plastic strains

    double* Epl = Ecg + 7;
    for (i = 0; i < 3; i++)
      Epl[i] += gam*nn[i];

    // Elasto-plastic tangent

    RealArray RHS(nvar*3,0.0);
    for (i = 0; i < 3; i++)
      RHS[nvar*i+i] = one2;
    if (!utl::solve(tres,RHS,&iPivot))
      return -4;
    for (a = 1; a <= 3; a++)
      for (b = a; b <= 3; b++)
        dtde(a,b) = RHS[nvar*(a-1)+b-1];
  }
  else // Elastic step (only tangent computation)
  {
    const RealArray& Ecg = be; // Elastic left Cauchy-Green tensor
    std::copy(Ecg.begin(),Ecg.end(),HVupd.begin());
    dtde(1,1) = dtde(2,2) = dtde(3,3) = one2*Bmod + two3*Smod;
    dtde(1,2) = dtde(1,3) = dtde(2,3) = one2*Bmod - one3*Smod;
  }


  // Compute the Cauchy stress

#if INT_DEBUG > 1
  std::cout <<"\nPlasticMaterial: dtde\n"<< dtde
            <<"\n                 tau = "<< tau << std::endl;
#endif

  RealArray Sig(6);
  for (i = 0; i < 6; i++)
    Sig[i] = (tau.x * nn_tr(p1[i],1)*nn_tr(p2[i],1) +
              tau.y * nn_tr(p1[i],2)*nn_tr(p2[i],2) +
              tau.z * nn_tr(p1[i],3)*nn_tr(p2[i],3)) / Jc;
  if (sigma.size() == 3)
    sigma = RealArray({ Sig[0], Sig[1], Sig[3] });
  else
    sigma = Sig;

  // Tangent transformation

  Matrix Cstm(6,6), Tmat(6,6);
  const Tensor& F = nn_tr;

  // Material tangent (computation in the principal basis)

  for (a = 1; a <= 3; a++)
  {
    // Upper 3x3 block of Cstm
    for (b = 1; b <= 3; b++)
      Cstm(b,a) = 2.0*dtde(b,a);
    Cstm(a,a)  -= 2.0*tau(a);

    // Lower 3x3 block of Cstm [ diagonal block ]

    b = a%3 + 1;
    double dtmp = ll2_tr(b) - ll2_tr(a);
    if (dtmp > tolb || dtmp < -tolb)
      Cstm(a+3,a+3) = (ll2_tr(a)*tau(b) - ll2_tr(b)*tau(a)) / dtmp;
    else
      Cstm(a+3,a+3) = dtde(a,a) - dtde(b,a) - tau(a);
  }

  /*
   * Form transformation matrix for a 4th rank tensor in matrix form
   *
   *    Tmat(a,b) = F(i,I)*F(j,J) : a -> I,J ; b -> i,j
   *
   *         a,b  |  1    2    3    4    5    6
   *        ------+-----------------------------
   *        (I,J) | 1,1  2,2  3,3  1,2  2,3  3,1
   *     or (i,j) |                2,1  3,2  1,3
   */

  auto i1 = [p1](int i) { return p1[i-1]; };
  auto i2 = [p2](int i) { return p2[i-1]; };

  for (i = 1; i <= 3; i++)
  {
    for (j = 1; j <= 3; j++)
      Tmat(i,j) = F(i1(j),i1(i)) * F(i2(j),i2(i));

    for (j = 4; j <= 6; j++)
      Tmat(i,j) = (F(i1(j),i1(i)) * F(i2(j),i2(i)) +
                   F(i2(j),i2(i)) * F(i1(j),i1(i))) * one2;
  }

  for (i = 4; i <= 6; i++)
  {
    for (j = 1; j <= 3; j++)
      Tmat(i,j) = (F(i1(j),i1(i)) * F(i2(j),i2(i)) +
                   F(i2(j),i2(i)) * F(i1(j),i1(i)));

    for (j = 4; j <= 6; j++)
      Tmat(i,j) = (F(i1(j),i1(i)) * F(i2(j),i2(i)) +
                   F(i2(j),i1(i)) * F(i1(j),i2(i)) +
                   F(i1(j),i2(i)) * F(i2(j),i1(i)) +
                   F(i2(j),i2(i)) * F(i1(j),i1(i))) * one2;
  }

#if INT_DEBUG > 1
  std::cout <<"PlasticMaterial::Cstm ="<< Cstm;
  std::cout <<"PlasticMaterial::Tmat ="<< Tmat;
#endif

  // Compute the triple matrix product: C = Tmat^t * Cstm * Tmat

  Matrix Ctmp;
  if (sigma.dim() == 3)
  {
    Ctmp.multiply(Cstm,Tmat);
    C.multiply(Tmat,Ctmp,true).multiply(1.0/Jc);
  }
  else // extract the 2D material matrix
  {
    C.multiply(Cstm,Tmat);
    Ctmp.multiply(Tmat,C,true).multiply(1.0/Jc);
    C.resize(sigma.size(),sigma.size());
    if (sigma.size() == 3)
      for (int i = 1; i <= 3; i++)
        for (int j = 1; j <= 3; j++)
          C(i,j) = Ctmp(i == 3 ? 4 : i, j == 3 ? 4 : j);
    else // 2D matrix with the zz-component included (plane strain)
      Ctmp.extractBlock(C,1,1);
  }

#if INT_DEBUG > 0
  std::cout <<"PlasticMaterial::sigma =\n"<< sigma;
  std::cout <<"PlasticMaterial::C ="<< C;
  std::cout <<"PlasticMaterial::HV(out) =";
  for (double v : HVc) std::cout <<" "<< v;
  std::cout << std::endl;
  std::cout.precision(oldPrec);
#endif

  const_cast<PlasticPoint*>(this)->updated = ierr > 0 ? 'd' : 'c';
  return ierr;
}


bool PlasticMaterial::PlasticPoint::yfunc (bool lIter, double Epp,
                                           double I1, double J2, double J3,
                                           double& f1, double& f2, double& f3,
                                           double& f11, double& f22,
                                           double& f33, double& f12,
                                           double& f13, double& f23,
                                           double& Ypr, double& YY,
                                           double& Yield) const
{
#if INT_DEBUG > 2
  std::cout <<"\nPlasticMaterial::yfunc: liter = "
            << std::boolalpha << lIter <<" Epp = "<< Epp
            <<"\n                        I1 = "<< I1 <<" J2 = "<< J2
            <<" J3 = "<< J3 << std::endl;
#endif
#ifndef epsZ
  const double epsZ  = 1.0e-16;
#endif
  const double one3  = 1.0/3.0;
  const double two3  = 2.0/3.0;
  const double sqt23 = sqrt(two3);

  f1 = f2 = f3 = f11 = f22 = f33 = f12 = f13 = f23 = 0.0;

  const int iYIELD  = static_cast<int>(pMAT[8]); // Type of yield function
  const double Hiso = sqt23 * pMAT[6];  // Isotropic hardening = sqrt(2/3)*H_iso
  const double Y0   = sqt23 * pMAT[9];  // Radius of yield = sqrt(2/3) Sig_y
  const double Yinf = sqt23 * pMAT[10]; // Radius of yield = sqrt(2/3) Sig_inf
  const double beta = pMAT[11];

  YY    = Y0 + Hiso*Epp;  // Sig_y(t_n) yield stress
  Ypr   = two3 * pMAT[6]; // Isotropic yield = 2/3 H_iso
  Yield = -YY;            // Default return (no yield)

  double aa, bb, c1, c2, q1, q2;
  double ss = J2 > 0.0 ? sqrt(2.0*J2) : 0.0;

  auto pow3 = [](double x){ return x*x*x; };

  switch (iYIELD) {

  case 1: // von Mises yield function (uses Y0, Yinf, H_iso, beta, J2)
    aa    = (Y0 - Yinf)*exp(-beta*Epp);
    YY    = Yinf + aa + Hiso*Epp;
    Ypr  -= sqt23*beta*aa;
    Yield = ss - YY;
    if (lIter && ss > epsZ)
    {
      // Compute yield function derivatives
      f2  = 1.0 / ss;
      f22 = -pow3(f2);
    }
    break;

  case 2: // Drucker-Prager yield function (uses Y0, Yinf, H_iso, I1, J2)
    aa    = Yinf;
    Yield = ss + one3*aa*I1 - YY;
    if (lIter)
    {
      // Compute yield function derivatives
      f1 = one3*aa;
      if (ss > epsZ)
      {
        f2  = 1.0 / ss;
        f22 = -pow3(f2);
      }
    }
    break;

  case 3: // Prager-Lode yield function (uses Y0, Yinf, H_iso, J2, J3)
    if (J2 < epsZ)
      Yield = -YY;
    else if (lIter)
    {
      double c1 = sqrt(0.5);
      double c2 = sqrt(13.5) * Yinf;

      Yield = ss + c2*J3/J2 - YY;

      // Compute yield function derivatives

      double J2_1  = 1.0 / J2;
      double J2_2  = J2_1*J2_1;
      double J2_12 = sqrt(J2_1);

      f2  =   c1*     J2_12      - c2*J3*     J2_2;
      f3  =                        c2   *     J2_1;
      f22 = - c1*pow3(J2_12)*0.5 + c2*J3*pow3(J2_1)*2.0;
      f23 =                      - c2   *     J2_2;
    }
    break;

  case 4: // Johnson-Cook yield function (uses Y0, Yinf, H_iso, J2)
    if (Epp < epsZ && -Epp < epsZ)
    {
      YY   = Y0;
      Ypr += sqt23*Yinf*pow(epsZ,beta-1.0);
    }
    else
    {
      YY   = Y0 + Yinf*pow(Epp,beta) + Hiso*Epp;
      Ypr += sqt23*beta*Yinf*pow(Epp,beta-1.0);
    }

    Yield = ss - YY;
    if (lIter && ss > epsZ)
    {
      // Compute yield function derivatives
      f2  =  1.0 / ss;
      f22 = -pow3(f2);
    }
    break;

  case 5: // Isotropic hardening curve (uses Hiso, J2)
    if (!hfn)
    {
      Yield = 0.0;
      std::cerr <<" *** PlasticMaterial::yfunc: No hardening function."
                << std::endl;
      return false;
    }

    YY    = sqt23 * hfn->eval(Epp);
    Ypr  += two3 * hfn->deriv(Epp);
    Yield = ss - YY;
    if (lIter && ss > epsZ)
    {
      // Compute yield function derivatives
      f2  =  1.0 / ss;
      f22 = -pow3(f2);
    }
    break;

  case 6: // Voce yield function (uses Y0, Hiso, c1, q1, c2, q2, J2)
    q1    = Yinf;
    c1    = beta;
    q2    = sqt23 * pMAT[13];
    c2    = pMAT[14];
    aa    = q1*exp(-c1*Epp);
    bb    = q2*exp(-c2*Epp);
    YY    = Y0 + q1 + q2 - aa - bb + Hiso*Epp;
    Ypr  += sqt23*(c1*aa + c2*bb);
    Yield = ss - YY;
    if (lIter && ss > epsZ)
    {
      // Compute yield function derivatives
      f2  =  1.0 / ss;
      f22 = -pow3(f2);
    }
    break;

  default:
    std::cerr <<" *** PlasticMaterial::yfunc: Unknown yield function "<< iYIELD
              << std::endl;
    return false;
  }

#if INT_DEBUG > 2
  std::cout <<"\nPlasticMaterial::yfunc: YY = "<< YY <<" Ypr = "<< Ypr
            <<" Yield = "<< Yield
            <<"\n                        F2 = "<< f2 <<" F22 = "<< f22
            << std::endl;
#endif
  return true;
}


double PlasticMaterial::PlasticPoint::energyIntegral (const SymmTensor& S,
                                                      const SymmTensor& E)
{
  Up += 0.5*ddot(S+Sp,E-Ep);
  Sp.copy(S);
  Ep.copy(E);
  return Up;
}


PlasticMaterial::ResultPoint::ResultPoint (const PlasticMaterial* prm,
                                           unsigned short int n)
  : PlasticPoint(prm,n)
{
  sigma_h = sigma_e = 0.0;
}


void PlasticMaterial::ResultPoint::principalStress (const SymmTensor& sigma)
{
  // Calculate some additional stress meassures
  sigma_h = sigma.trace() / double(sigma.size() > 3 ? 3 : sigma.dim());
  sigma_e = sigma.vonMises();
  sigma.principal(prin);
}


double PlasticMaterial::ResultPoint::getTriax () const
{
  return sigma_e > 0.0 ? sigma_h / sigma_e : 0.0;
}


double PlasticMaterial::ResultPoint::getLode () const
{
  double denom = prin.x - prin.z;
  return denom > 0.0 ? (2.0*prin.y - prin.x - prin.z) / denom : 0.0;
}
