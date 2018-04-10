// $Id$
//==============================================================================
//!
//! \file KirchhoffLoveShell.C
//!
//! \date Feb 25 2018
//!
//! \author ... and ... / NTNU
//!
//! \brief Class for linear Kirchhoff-Love thin shell problems.
//!
//==============================================================================

#include "KirchhoffLoveShell.h"
#include "LinIsotropic.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Vec3Oper.h"
#include "Utilities.h"
#include "IFEM.h"


KirchhoffLoveShell::KirchhoffLoveShell () : KirchhoffLove(3)
{
  tracFld = nullptr;
  fluxFld = nullptr;
}


void KirchhoffLoveShell::printLog () const
{
  IFEM::cout <<"KirchhoffLoveShell: thickness = "<< thickness
             <<", gravity = "<< gravity << std::endl;

  if (!material)
  {
    static LinIsotropic defaultMat;
    const_cast<KirchhoffLoveShell*>(this)->material = &defaultMat;
  }

  material->printLog();
}


bool KirchhoffLoveShell::evalInt (LocalIntegral& elmInt,
                                  const FiniteElement& fe,
                                  const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (eM) // Integrate the mass matrix
    this->formMassMatrix(elMat.A[eM-1],fe.N,X,fe.detJxW);

  if (eK) // Integrate the stiffness matrix
    this->evalK(elMat.A[eK-1],fe,X);

  if (eS) // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(elMat.b[eS-1],fe.N,fe.iGP,X,fe.detJxW);

  return true;
}


void KirchhoffLoveShell::formMassMatrix (Matrix& EM, const Vector& N,
                                         const Vec3& X, double detJW) const
{
  double rhow = material->getMassDensity(X)*thickness*detJW;
  if (rhow == 0.0) return;

  for (size_t a = 1; a <= N.size(); a++)
    for (size_t b = 1; b <= N.size(); b++)
      for (unsigned short int i = 1; i <= 3; i++)
        EM(3*(a-1)+i,3*(b-1)+i) += rhow*N(a)*N(b);
}


void KirchhoffLoveShell::formBodyForce (Vector& ES, const Vector& N, size_t iP,
                                        const Vec3& X, double detJW) const
{
  double p = this->getPressure(X);
  if (p == 0.0) return;

  for (size_t a = 1; a <= N.size(); a++)
    ES(3*a) += N(a)*p*detJW;

  // Store pressure value for visualization
  if (iP < presVal.size())
    presVal[iP] = std::make_pair(X,Vec3(0.0,0.0,p));
}


Vec3 KirchhoffLoveShell::getTraction (const Vec3& X, const Vec3& n) const
{
  if (fluxFld)
    return (*fluxFld)(X);
  else if (tracFld)
    return (*tracFld)(X,n);
  else
    return Vec3();
}


bool KirchhoffLoveShell::evalK (Matrix& EK, const FiniteElement& fe,
                                const Vec3& X) const
{

    int ndof = fe.N.size()*3;

      Matrix kem(ndof,ndof);
      Matrix keb(ndof,ndof);

      Matrix dN_ca(3,ndof);
      Matrix dM_ca(3,ndof);

      Matrix Dm;
      Matrix Db;
      if (!this->formDmatrix(Dm,Db,fe,X))
        return false;

      Matrix Bm;
      Matrix Bb;
      if (!this->formBmatrix(Bm,Bb,fe))
          return false;

      dN_ca.multiply(Dm,Bm);
      dM_ca.multiply(Db,Bb);

      kem.multiply(Bm,dN_ca,true,false,true);
      keb.multiply(Bb,dM_ca,true,false,true);
      kem.multiply(fe.detJxW);
      keb.multiply(fe.detJxW);

      EK.add(kem).add(keb);


      return true;

}

bool KirchhoffLoveShell::formDmatrix (Matrix& Dm, Matrix& Db, const FiniteElement& fe, const Vec3& X, bool invers) const
{
    SymmTensor dummy(2); double U;
    Matrix D;
    if (!material->evaluate(D,dummy,U,fe,X,dummy,dummy, invers ? -1 : 1))
        return false;

    Matrix Dtemp = D;
    double factorm = thickness;
    double factorb = thickness*thickness*thickness/12;

    Dm = D.multiply(invers ? 1.0/factorm : factorm);
    Db = Dtemp.multiply(invers ? 1.0/factorb : factorb);

    return true;
}

bool KirchhoffLoveShell::formBmatrix (Matrix& Bm, Matrix& Bb, const FiniteElement& fe) const
{
    // Covariant basis vectors
    Vec3 g1 = fe.G.getColumn(1);
    Vec3 g2 = fe.G.getColumn(2);

    // Basis vector g3
    Vec3 g3(g1,g2);
    double lg3 = g3.length();
    Vec3 n(g3);
    n.normalize();

    // Covariant metric gab
    double gab11 = g1*g1;
    double gab12 = g1*g2;
    double gab22 = g2*g2;

    // Contravariant metric gab_con and base vectors g_con
    double invdetgab =  1.0/(gab11*gab22-gab12*gab12);
    double gab_con11 =  invdetgab*gab22;
    double gab_con12 = -invdetgab*gab12;
    double gab_con22 =  invdetgab*gab11;
    Vec3 g1_con = g1*gab_con11 + g2*gab_con12;
    Vec3 g2_con = g1*gab_con12 + g2*gab_con22;

    // Local cartesian coordinates
    Vec3 e1(g1);     e1.normalize();
    Vec3 e2(g2_con); e2.normalize();

    // Transformation matrix T from contravariant to local cartesian basis
    Matrix T(3,3);
    double eg11 = e1*g1_con;
    double eg12 = e1*g2_con;
    double eg21 = e2*g1_con;
    double eg22 = e2*g2_con;
    T(1,1) =      eg11*eg11;
    T(1,2) =      eg12*eg12;
    T(1,3) = 2.0* eg11*eg12;
    T(2,1) =      eg21*eg21;
    T(2,2) =      eg22*eg22;
    T(2,3) = 2.0* eg21*eg22;
    T(3,1) = 2.0* eg11*eg21;
    T(3,2) = 2.0* eg12*eg22;
    T(3,3) = 2.0*(eg11*eg22+eg12*eg21);

    // Strain
    int ndof = fe.N.size()*3;
    Matrix dE_ca(3,ndof); // dE_ca = epsilon cartesian coordinates
    Matrix dK_ca(3,ndof); // dK_ca = kappa cartesian coordinates
    Matrix dE_cu(3,ndof); // dE_cu = epsilon curvelinear coordinate system
    Matrix dK_cu(3,ndof); // dK_cu = kappa curvelinear

      for (int i = 1; i <= ndof; i++)
      {
        double dummy = i;
        double k = ceil(dummy/3);
        int dir = i-3*(k-1);

        dE_cu(1,i) = fe.dNdX(k,1)*g1(dir);
        dE_cu(2,i) = fe.dNdX(k,2)*g2(dir);
        dE_cu(3,i) = 0.5*(fe.dNdX(k,1)*g2(dir) + fe.dNdX(k,2)*g1(dir));

        Vec3 dg1;
        dg1(dir) = fe.dNdX(k,1);

        Vec3 dg2;
        dg2(dir) = fe.dNdX(k,2);

        Vec3 dummy2(dg1,g2);
        Vec3 dg3(g1,dg2);
        dg3 = dg3 + dummy2;

        double g3dg3lg3_3 = g3*dg3/(lg3*lg3*lg3); // eq 5.31 last part

        Vec3 dn; // eq 5.31 dn = a3,rs (kanskje)
        for (int kaffi = 1; kaffi <= 3; kaffi ++)
          {
             dn(kaffi) = dg3(kaffi)/lg3 - g3(kaffi)*g3dg3lg3_3;
          }



        dK_cu(1,i) = -(fe.d2NdX2(k,1,1)*n(dir) + fe.H.getColumn(1)*dn);
        dK_cu(2,i) = -(fe.d2NdX2(k,2,2)*n(dir) + fe.H.getColumn(2)*dn);
        dK_cu(3,i) = -(fe.d2NdX2(k,1,2)*n(dir) + fe.H.getColumn(3)*dn);


      } // for int i = 0; i <= ndof; i++

      Bm.multiply(T,dE_cu); // Bm = dE_ca.multiply(T,dE_cu);
      Bb.multiply(T,dK_cu); // Bb = dK_ca.multiply(T,dK_cu);
    return true;

}


bool KirchhoffLoveShell::evalBou (LocalIntegral& elmInt,
                                  const FiniteElement& fe,
                                  const Vec3& X, const Vec3& normal) const
{
  if (!eS)
  {
    std::cerr <<" *** KirchhoffLoveShell::evalBou: No load vector."<< std::endl;
    return false;
  }
  else if (!fluxFld && !tracFld)
  {
    std::cerr <<" *** KirchhoffLoveShell::evalBou: No tractions."<< std::endl;
    return false;
  }

  Vec3 T = this->getTraction(X,normal);
  Vector& ES = static_cast<ElmMats&>(elmInt).b[eS-1];
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (unsigned short int i = 1; i <= 3; i++)
      ES(3*(a-1)+i) += T[i-1]*fe.N(a)*fe.detJxW;

  return true;
}


bool KirchhoffLoveShell::evalSol (Vector& s,
                                  const FiniteElement& fe, const Vec3& X,
                                  const std::vector<int>& MNPC) const
{
  // Extract element displacements
  Vector eV;
  if (!primsol.empty() && !primsol.front().empty())
  {
    int ierr = utl::gather(MNPC,3,primsol.front(),eV);
    if (ierr > 0)
    {
      std::cerr <<" *** KirchhoffLoveShell::evalSol: Detected "
                << ierr <<" node numbers out of range."<< std::endl;
      return false;
    }
  }
  s.reserve(18);

  // Evaluate the stress resultant tensor
  Vector sb;
  if (!this->evalSol(s,sb,eV,fe,X,true))
    return false;

  // Append the bending moment to the stress array
  s.insert(s.end(),sb.begin(),sb.end());

  // Calculate top and bottom surface stresses
  SymmTensor sigma(2);
  Vec3       sigma_p;
  for (int isurf = -1; isurf < 2; isurf += 2)
  {
    double* p = const_cast<double*>(sigma.ptr());
    for (size_t i = 0; i < 3; i++)
      p[i] = (s[i] - isurf*sb[i]*6.0/thickness)/thickness;

    sigma.principal(sigma_p);
    s.insert(s.end(),sigma.ptr(),sigma.ptr()+3);
    s.insert(s.end(),sigma_p.ptr(),sigma_p.ptr()+2);
    s.insert(s.end(),sigma.vonMises());
  }

  return true;
}


bool KirchhoffLoveShell::evalSol (Vector& sm, Vector& sb, const Vector& eV,
                                  const FiniteElement& fe, const Vec3& X,
                                  bool toLocal) const
{
  if (eV.empty())
  {
    std::cerr <<" *** KirchhoffLoveShell::evalSol: No displacement vector."
              << std::endl;
    return false;
  }
  else if (eV.size() != 3*fe.d2NdX2.dim(1))
  {
    std::cerr <<" *** KirchhoffLoveShell::evalSol: Invalid displacement vector."
              <<"\n     size(eV) = "<< eV.size() <<"   size(d2NdX2) = "
              << fe.d2NdX2.dim(1) <<","<< fe.d2NdX2.dim(2)*fe.d2NdX2.dim(3)
              << std::endl;
    return false;
  }

  // Compute the strain-displacement matrices Bm and Bb
  Matrix Bm, Bb;
  if (!this->formBmatrix(Bm,Bb,fe))
    return false;

  // Evaluate the constitutive matrices at this point
  Matrix Dm, Db;
  if (!this->formDmatrix(Dm,Db,fe,X))
    return false;

  // Evaluate the membran strain and curvature tensors
  SymmTensor epsilon(2), kappa(2);
  if (!Bm.multiply(eV,epsilon)) // epsilon = B*eV
    return false;
  if (!Bb.multiply(eV,kappa)) // kappa = B*eV
    return false;

  // Evaluate the stress resultant tensors
  SymmTensor n(2), m(2);
  if (!Dm.multiply(epsilon,n)) // n = Dm*epsilon
    return false;
  if (!Db.multiply(-1.0*kappa,m)) // m = -Db*kappa
    return false;

  // Congruence transformation to local coordinate system at current point
  if (toLocal && locSys)
  {
    n.transform(locSys->getTmat(X));
    m.transform(locSys->getTmat(X));
  }

  sm = n;
  sb = m;
  return true;
}


std::string KirchhoffLoveShell::getField1Name (size_t i,
                                               const char* prefix) const
{
  if (i >= 3) return "";

  char name = 'u'+i;
  if (!prefix)
    return std::string(1,name);

  return prefix + std::string(" ") + std::string(1,name);
}


std::string KirchhoffLoveShell::getField2Name (size_t i,
                                               const char* prefix) const
{
  static const char* s[12] = { "n_xx", "n_yy", "n_xy", "m_xx", "m_yy", "m_xy",
                               "sigma_x", "sigma_y", "tau_xy",
                               "sigma_1", "sigma_2", "sigma_m" };

  std::string name(s[i < 6 ? i : 6 + i%6]);
  if (i >= 12)
    name = "Top " + name;
  else if (i >= 6)
    name = "Bottom " + name;

  if (!prefix)
    return name;

  return prefix + std::string(" ") + name;
}


NormBase* KirchhoffLoveShell::getNormIntegrand (AnaSol*) const
{
  return new KirchhoffLoveShellNorm(*const_cast<KirchhoffLoveShell*>(this));
}


KirchhoffLoveShellNorm::KirchhoffLoveShellNorm (KirchhoffLoveShell& p)
  : NormBase(p)
{
  nrcmp = 6;
}


bool KirchhoffLoveShellNorm::evalInt (LocalIntegral& elmInt,
                                      const FiniteElement& fe,
                                      const Vec3& X) const
{
  KirchhoffLoveShell& problem = static_cast<KirchhoffLoveShell&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the inverse constitutive matrices at this point
  Matrix Dm, Db;
  if (!problem.formDmatrix(Dm,Db,fe,X,true))
    return false;

  // Evaluate the finite element stress field
  Vector mh, nh, errm, errn;
  if (!problem.evalSol(nh,mh,pnorm.vec.front(),fe,X))
    return false;

  size_t ip = 0;

  // Integrate the energy norm a(u^h,u^h)
  pnorm[ip++] += (nh.dot(Dm*nh) + mh.dot(Db*mh))*fe.detJxW;

  // Evaluate the body load
  double p = problem.getPressure(X);
  // Evaluate the transverse displacement field
  double w = pnorm.vec.front().dot(fe.N,2,3);
  // Integrate the external energy (p,u^h)
  pnorm[ip++] += p*w*fe.detJxW;

  // Integrate the area
  pnorm[ip++] += fe.detJxW;

#if INT_DEBUG > 3
  std::cout <<"KirchhoffLovePlateNorm::evalInt("<< fe.iel <<", "<< X <<"):";
#endif

  for (const Vector& psol : pnorm.psol)
    if (!psol.empty())
    {
      // Evaluate the projected solution
      Vector nr(3), mr(3);
      for (unsigned short int j = 0; j < 6; j++)
        if (j < 3)
          nr[j] = psol.dot(fe.N,j,nrcmp);
        else
          mr[j-3] = psol.dot(fe.N,j,nrcmp);

#if INT_DEBUG > 3
      std::cout <<"\n\ts^r =";
      for (double v : nr) std::cout <<" "<< v;
      for (double v : mr) std::cout <<" "<< v;
#endif

      // Integrate the energy norm a(u^r,u^r)
      pnorm[ip++] += (nr.dot(Dm*nr) + mr.dot(Db*mr))*fe.detJxW;

      // Integrate the error in energy norm a(u^r-u^h,u^r-u^h)
      errn = nr - nh;
      errm = mr - mh;
      pnorm[ip++] += (errn.dot(Dm*errn) + errm.dot(Db*errm))*fe.detJxW;

      // Integrate the L2-norms (n^r,n^r) and (m^r,m^r)
      pnorm[ip++] += nr.dot(nr)*fe.detJxW;
      pnorm[ip++] += mr.dot(mr)*fe.detJxW;
      // Integrate the error in L2-norm (m^r-m^h,m^r-m^h)
      pnorm[ip++] += errn.dot(errn)*fe.detJxW;
      pnorm[ip++] += errm.dot(errm)*fe.detJxW;
    }

  if (ip == pnorm.size())
    return true;

  std::cerr <<" *** KirchhoffLoveShellNorm::evalInt: Internal error, ip="
            << ip <<" != pnorm.size()="<< pnorm.size() << std::endl;
  return false;
}


bool KirchhoffLoveShellNorm::evalBou (LocalIntegral& elmInt,
                                      const FiniteElement& fe,
                                      const Vec3& X, const Vec3& normal) const
{
  KirchhoffLoveShell& problem = static_cast<KirchhoffLoveShell&>(myProblem);
  if (!problem.haveLoads()) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the surface traction and displacement field
  Vec3 u, T = problem.getTraction(X,normal);
  for (int i = 0; i < 3; i++)
    u[i] = pnorm.vec.front().dot(fe.N,i,3);

  // Integrate the external energy
  pnorm[1] += T*u*fe.detJxW;
  return true;
}


int KirchhoffLoveShellNorm::getIntegrandType () const
{
  return SECOND_DERIVATIVES;
}


size_t KirchhoffLoveShellNorm::getNoFields (int group) const
{
  if (group == 0)
    return this->NormBase::getNoFields();
  else if (group == 1 || group == -1)
    return 3;
  else if (group > 0 || !prjsol[-group-2].empty())
    return 6;
  else
    return 0;
}


std::string KirchhoffLoveShellNorm::getName (size_t i, size_t j,
                                             const char* prefix) const
{
  if (i == 0 || j == 0 || j > 6 || (i == 1 && j > 3))
    return this->NormBase::getName(i,j,prefix);

  static const char* u[3] = {
    "a(u^h,u^h)^0.5",
    "(p,u^h)^0.5",
    "area"
  };

  static const char* p[6] = {
    "a(u^r,u^r)^0.5",
    "a(e,e)^0.5, e=u^r-u^h",
    "(n^r,n^r)^0.5",
    "(m^r,m^r)^0.5",
    "(e,e)^0.5, e=n^r-n^h",
    "(e,e)^0.5, e=m^r-m^h",
  };

  std::string name(i > 1 ? p[j-1] : u[j-1]);
  return prefix ? prefix + std::string(" ") + name : name;
}
