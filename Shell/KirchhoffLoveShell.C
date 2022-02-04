// $Id$
//==============================================================================
//!
//! \file KirchhoffLoveShell.C
//!
//! \date Feb 25 2018
//!
//! \author Simen Skogholt Haave and Marit Gaarder Rakvaag / NTNU
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
#include "IFEM.h"


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
    this->formBodyForce(elMat.b[eS-1],fe.N,fe.iGP,X,
                        this->getShellNormal(fe.G),fe.detJxW);

  if (gS && !presFld.empty()) // Integrate the pressure load vector gradient
    this->formBodyForce(elMat.b[gS-1],fe.N,fe.iGP,X,
                        this->getShellNormal(fe.G),fe.detJxW,true);

  return true;
}


bool KirchhoffLoveShell::evalK (Matrix& EK, const FiniteElement& fe,
                                const Vec3& X) const
{
  // Constitutive matrices
  Matrix Dm, Db;
  if (!this->formDmatrix(Dm,Db,fe,X))
    return false;

  // Strain-displacement matrices
  Matrix Bm, Bb;
  if (!this->formBmatrix(Bm,Bb,fe))
    return false;

  Matrix dN_ca, dM_ca;
  dN_ca.multiply(Dm,Bm); // dN_ca = Dm*Bm, cartesian coordinates
  dM_ca.multiply(Db,Bb); // dM_ca = Db*Bb, cartesian coordinates

  // Element stiffness matrices,
  // with membrane and bending contributions separated
  Matrix kem, keb;
  kem.multiply(Bm,dN_ca,true); // kem = Bm^t*dN_ca = Bm^t*Dm*Bm
  keb.multiply(Bb,dM_ca,true); // keb = Bb^t*dM_ca = Bb^t*Db*Bb

  EK.add(kem,fe.detJxW).add(keb,fe.detJxW); // EK += (kem+keb)*|J|*w

  return true;
}


bool KirchhoffLoveShell::formDmatrix (Matrix& Dm, Matrix& Db,
                                      const FiniteElement& fe, const Vec3& X,
                                      bool invers) const
{
  SymmTensor dummy(2); double U;
  if (!material->evaluate(Dm,dummy,U,fe,X,dummy,dummy, invers ? -1 : 1))
    return false;

  double factorm = thickness;
  double factorb = thickness*thickness*thickness/12.0;

  Db = Dm;
  Dm.multiply(invers ? 1.0/factorm : factorm);
  Db.multiply(invers ? 1.0/factorb : factorb);

  return true;
}


bool KirchhoffLoveShell::formBmatrix (Matrix& Bm, Matrix& Bb,
                                      const FiniteElement& fe) const
{
  // Calculate metrics
  Vec3 g1, g2, g3, n, gab; Matrix T(3,3);
  double lg3 = this->getMetrics(fe.G,g1,g2,g3,n,gab,&T);
  double lg3_2 = lg3*lg3;
  size_t nenod = fe.dNdX.rows();
  size_t nedof = 3*nenod;

  // Strains
  Matrix dE_cu(3,nedof); // dE_cu = epsilon curvelinear coordinate system
  Matrix dK_cu(3,nedof); // dK_cu = kappa curvelinear

  for (size_t k = 1; k <= nenod; k++)
    for (int dir = 1; dir <= 3; dir++)
    {
      Vec3 dg1, dg2;
      dg1(dir) = fe.dNdX(k,1);
      dg2(dir) = fe.dNdX(k,2);
      Vec3 dg3 = Vec3(g1,dg2) + Vec3(dg1,g2);
      Vec3 dn  = (dg3 - g3*((g3*dg3)/lg3_2))/lg3;
      Vec3 bv  = dn * fe.H;
      size_t i = 3*k-3 + dir;

      dE_cu(1,i) = fe.dNdX(k,1)*g1(dir);
      dE_cu(2,i) = fe.dNdX(k,2)*g2(dir);
      dE_cu(3,i) = 0.5*(fe.dNdX(k,1)*g2(dir) + fe.dNdX(k,2)*g1(dir));

      dK_cu(1,i) = -fe.d2NdX2(k,1,1)*n(dir) - bv.x;
      dK_cu(2,i) = -fe.d2NdX2(k,2,2)*n(dir) - bv.y;
      dK_cu(3,i) = -fe.d2NdX2(k,1,2)*n(dir) - bv.z;
    }

  Bm.multiply(T,dE_cu); // Bm = dE_ca = T * dE_cu
  Bb.multiply(T,dK_cu); // Bb = dK_ca = T * dK_cu

  return true;
}


double KirchhoffLoveShell::getMetrics (const Matrix& G,
                                       Vec3& g1, Vec3& g2, Vec3& g3, Vec3& n,
                                       Vec3& gab, Matrix* Tlc) const
{
  // Covariant basis vectors
  g1 = G.getColumn(1);
  g2 = G.getColumn(2);
  n = g3.cross(g1,g2);

  // Covariant metric gab
  gab.x = g1*g1;
  gab.y = g2*g2;
  gab.z = g1*g2;

  if (Tlc)
  {
    Matrix& T = *Tlc;

    // Contravariant metric gab_con and base vectors g_con
    double invdetgab =  1.0/(gab.x*gab.y - gab.z*gab.z);
    double gab_con11 =  invdetgab*gab.y;
    double gab_con12 = -invdetgab*gab.z;
    double gab_con22 =  invdetgab*gab.x;
    Vec3 g1_con = g1*gab_con11 + g2*gab_con12;
    Vec3 g2_con = g1*gab_con12 + g2*gab_con22;

    // Local cartesian coordinates
    Vec3 e1(g1);     e1.normalize();
    Vec3 e2(g2_con); e2.normalize();

    // Transformation matrix from contravariant to local cartesian basis
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
  }

  return n.normalize();
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
    for (int i = 1; i <= 3; i++)
      ES(3*(a-1)+i) += T[i-1]*fe.N(a)*fe.detJxW;

  // Store traction value for visualization
  if (fe.iGP < tracVal.size() && !T.isZero())
  {
    tracVal[fe.iGP].first = X;
    tracVal[fe.iGP].second += T;
  }

  if (gS)
  {
    T = this->getTraction(X,normal,true);
    Vector& GS = static_cast<ElmMats&>(elmInt).b[gS-1];
    for (size_t a = 1; a <= fe.N.size(); a++)
      for (int i = 1; i <= 3; i++)
        GS(3*(a-1)+i) += T[i-1]*fe.N(a)*fe.detJxW;
  }

  return true;
}


bool KirchhoffLoveShell::evalSol (Vector& s, const Vectors& eV,
                                  const FiniteElement& fe, const Vec3& X,
                                  bool toLocal) const
{
  // Evaluate the stress resultants (in-plane forces and bending moments)
  s.reserve(18); Vector sb;
  if (!this->evalSol(s,sb,eV,fe,X,toLocal))
    return false;
  else
    s.insert(s.end(),sb.begin(),sb.end());

  // Calculate top and bottom surface stresses;
  // stress tensor components, principal stresses and von Mises stress
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


bool KirchhoffLoveShell::evalSol (Vector& sm, Vector& sb, const Vectors& eV,
                                  const FiniteElement& fe, const Vec3& X,
                                  bool toLocal) const
{
  if (eV.empty() || eV.front().empty())
  {
    std::cerr <<" *** KirchhoffLoveShell::evalSol: No displacement vector."
              << std::endl;
    return false;
  }
  else if (eV.front().size() != 3*fe.d2NdX2.dim(1))
  {
    std::cerr <<" *** KirchhoffLoveShell::evalSol: Invalid displacement vector."
              <<"\n     size(eV) = "<< eV.front().size() <<"   size(d2NdX2) = "
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

  // Evaluate the membrane strain and curvature tensors
  SymmTensor epsilon(2), kappa(2);
  if (!Bm.multiply(eV.front(),epsilon)) // epsilon = B*eV
    return false;
  if (!Bb.multiply(eV.front(),kappa)) // kappa = B*eV
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


Vec3 KirchhoffLoveShell::getShellNormal (const Matrix& G) const
{
  Vec3 n(G.getColumn(1),G.getColumn(2));
  n.normalize();
  return n;
}


void KirchhoffLoveShell::primaryScalarFields (Matrix& field)
{
  if (field.rows() != 3) return;

  // Insert the absolute value as the fourth solution component
  double* u = new double[field.cols()];
  for (size_t c = 1; c <= field.cols(); c++)
    u[c-1] = field.getColumn(c).norm2();
  field.expandRows(1);
  field.fillRow(4,u);
  delete[] u;
}


std::string KirchhoffLoveShell::getField1Name (size_t i,
                                               const char* prefix) const
{
  if (i == 3)
  {
    std::string name("sqrt(u^2+v^2+w^2)");
    if (prefix)
      name = prefix + std::string(" ") + name;
    return name;
  }
  else if (i > 3)
  {
    if (prefix)
      return prefix + std::string(" displacement");
    return "displacement";
  }

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
  nrcmp = myProblem.getNoFields(2);
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
  if (!problem.evalSol(nh,mh,pnorm.vec,fe,X))
    return false;

  // Evaluate the pressure load and displacement field
  Vec3 u, p = problem.getPressure(X,problem.getShellNormal(fe.G));
  for (int i = 0; i < 3; i++)
    u[i] = pnorm.vec.front().dot(fe.N,i,3);

  // Integrate the energy norm a(u^h,u^h)
  pnorm[0] += (nh.dot(Dm*nh) + mh.dot(Db*mh))*fe.detJxW;
  // Integrate the external energy (p,u^h)
  pnorm[1] += p*u*fe.detJxW;

  if (pnorm.psol.empty() && pnorm.size() == 2)
    return true; // no projection in this run

#if INT_DEBUG > 3
  std::cout <<"KirchhoffLovePlateNorm::evalInt("<< fe.iel <<", "<< X;
  std::cout <<"):\n\ts^h =";
  for (double v : nh) std::cout <<" "<< v;
  for (double v : mh) std::cout <<" "<< v;
#endif

  size_t ip = 2;
  for (const Vector& psol : pnorm.psol)
    if (!psol.empty())
    {
      // Evaluate the projected solution
      Vector nr(3), mr(3);
      for (size_t j = 0; j < 6; j++)
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

#if INT_DEBUG > 3
  std::cout << std::endl;
#endif
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
  const KirchhoffLoveShell& problem = static_cast<const KirchhoffLoveShell&>(myProblem);
  if (!problem.haveLoads('B')) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the surface traction and displacement field
  Vec3 u, T = problem.getTraction(X,normal);
  for (int i = 0; i < 3; i++)
    u[i] = pnorm.vec.front().dot(fe.N,i,3);

  // Integrate the external energy
  pnorm[1] += T*u*fe.detJxW;
  return true;
}


size_t KirchhoffLoveShellNorm::getNoFields (int group) const
{
  if (group == 0)
    return this->NormBase::getNoFields();
  else if (group == 1 || group == -1)
    return 2;
  else if (group > 0 || !prjsol[-group-2].empty())
    return 6;
  else
    return 0;
}


std::string KirchhoffLoveShellNorm::getName (size_t i, size_t j,
                                             const char* prefix) const
{
  if (i == 0 || j == 0 || j > 6 || (i == 1 && j > 2))
    return this->NormBase::getName(i,j,prefix);

  static const char* u[2] = {
    "a(u^h,u^h)^0.5",
    "(p,u^h)^0.5"
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
