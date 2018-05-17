// $Id$
//==============================================================================
//!
//! \file NLKirchhoffLoveShell.C
//!
//! \date Apr 20 2018
//!
//! \author Simen Skogholt Haave and Marit Gaarder Rakvaag / NTNU
//!
//! \brief Class for nonlinear Kirchhoff-Love thin shell problems.
//!
//==============================================================================

#include "NLKirchhoffLoveShell.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "CoordinateMapping.h"
#include "Vec3Oper.h"


void NLKirchhoffLoveShell::setMode (SIM::SolutionMode mode)
{
  this->KirchhoffLoveShell::setMode(mode);
  if (mode == SIM::STATIC || mode == SIM::RHS_ONLY)
  {
    iS = 1;
    primsol.resize(1);
  }
}


int NLKirchhoffLoveShell::getIntegrandType () const
{
  return SECOND_DERIVATIVES | UPDATED_NODES;
}


bool NLKirchhoffLoveShell::evalInt (LocalIntegral& elmInt,
                                    const FiniteElement& fe,
                                    const Vec3& X) const
{
  Matrix Gd, Hd;
  if (elmInt.vec.size() > 1)
  {
    // Co-variant basis and Hessian in deformed configuration
    Matrix3D Hess;
    Gd.multiplyMat(elmInt.vec.back(),fe.dNdX);
    if (Hess.multiplyMat(elmInt.vec.back(),fe.d2NdX2))
      utl::Hessian(Hess,Hd);
    else
      return false;
  }

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (eM) // Integrate the mass matrix
    this->formMassMatrix(elMat.A[eM-1],fe.N,X,fe.detJxW);

  if (eK && iS) // Integrate the stiffness matrix and internal forces
    this->evalKandS(elMat.A[eK-1],elMat.b[iS-1],fe,fe.G,Gd,fe.H,Hd,X);
  else if (iS) // Integrate the internal forces only
  {
    Matrix dummyEK;
    this->evalKandS(dummyEK,elMat.b[iS-1],fe,fe.G,Gd,fe.H,Hd,X);
  }

  if (eS) // Integrate the load vector due to gravitation and other body forces
  {
    Vec3 n;
    if (presFld)
      n = this->getShellNormal(Gd.empty() ? fe.G : Gd);
    this->formBodyForce(elMat.b[eS-1],fe.N,fe.iGP,X,n,fe.detJxW);
  }

  return true;
}


bool NLKirchhoffLoveShell::evalKandS (Matrix& EK, Vector& ES,
                                      const FiniteElement& fe,
                                      const Matrix& G0, const Matrix& Gn,
                                      const Matrix& H0, const Matrix& Hn,
                                      const Vec3& X) const
{
  Matrix Dm, Db;
  if (!this->formDmatrix(Dm,Db,fe,X))
    return false;

  // Calculate the metrics
  Vec3 g1, g2, g3, n, gab, Gab; Matrix T(3,3);
  double lg3 = this->getMetrics(G0,g1,g2,g3,n,Gab,&T);
  Vec3 bv, Bv = n * H0;

  if (Gn.empty())
  {
    // Initial configuration, no deformation yet
    gab = Gab;
    bv  = Bv;
  }
  else
  {
    // Deformed configuration
    lg3 = this->getMetrics(Gn,g1,g2,g3,n,gab);
    bv  = n * Hn;
  }

#if INT_DEBUG > 1
  std::cout <<"\nNLKirchhoffLoveShell::evalKandS(X="<< X <<", iGP="<< fe.iGP
            <<")\n\tg1 = "<< g1 <<"\n\tg2 = "<< g2 <<"\n\tg3 = "<< g3
            <<"\n\tn = "<< n <<"\n\tgab = "<< gab <<"\n\tgab = "<< bv <<"\n";
#endif

  // Strain vectors referred to curvilinear coordinate system
  Vec3 E_cu = 0.5*(gab-Gab);
  Vec3 K_cu = (Bv - bv);

  // Strain vectors referred to cartesian coordinate system
  Vec3 E_ca = T * E_cu;
  Vec3 K_ca = T * K_cu;

  // Stress resultants referred to cartesian coordinate system
  Vector N_ca, M_ca;
  if (!Dm.multiply(E_ca.vec(),N_ca)) // N_ca = Dm*E_ca
    return false;
  if (!Db.multiply(K_ca.vec(),M_ca)) // M_ca = Db*K_ca
    return false;

#if INT_DEBUG > 1
  std::cout <<"\tE_ca = "<< E_ca <<"\n\tK_ca = "<< K_ca
            <<"\nN_ca:"<< N_ca <<"M_ca:"<< M_ca;
#endif

  // Establish the strain-displacement matrices
  size_t nenod = fe.dNdX.rows();
  size_t nedof = 3*nenod;
  Matrix dE_cu(3,nedof); // epsilon curvelinear coordinate system
  Matrix dK_cu(3,nedof); // kappa curvelinear coordinate system
  Matrix dg3(3,nedof), dn(3,nedof);
  Vector g3dg3lg3_3(nedof), g3dg3(nedof);
  double lg3_3 = lg3*lg3*lg3;
  double lg3_5 = lg3_3*lg3*lg3;
  for (size_t k = 1; k <= nenod; k++)
    for (int dir = 1; dir <= 3; dir++)
    {
      int i = 3*k-3 + dir;

      dE_cu(1,i) = fe.dNdX(k,1)*g1(dir);
      dE_cu(2,i) = fe.dNdX(k,2)*g2(dir);
      dE_cu(3,i) = 0.5*(fe.dNdX(k,1)*g2(dir) + fe.dNdX(k,2)*g1(dir));

      Vec3 dg1, dg2;
      dg1(dir) = fe.dNdX(k,1);
      dg2(dir) = fe.dNdX(k,2);
      dg3(1,i) = dg1(2)*g2(3) - dg1(3)*g2(2) + g1(2)*dg2(3) - g1(3)*dg2(2);
      dg3(2,i) = dg1(3)*g2(1) - dg1(1)*g2(3) + g1(3)*dg2(1) - g1(1)*dg2(3);
      dg3(3,i) = dg1(1)*g2(2) - dg1(2)*g2(1) + g1(1)*dg2(2) - g1(2)*dg2(1);
      g3dg3(i) = g3*dg3.getColumn(i);
      g3dg3lg3_3(i) = g3dg3(i)/lg3_3;
      Vec3 dni = dg3.getColumn(i)/lg3 - g3*g3dg3lg3_3(i);
      Vec3 dbv = dni * (Hn.empty() ? H0 : Hn);
      dn.fillColumn(i,dni.vec());

      dK_cu(1,i) = -fe.d2NdX2(k,1,1)*n(dir) - dbv.x;
      dK_cu(2,i) = -fe.d2NdX2(k,2,2)*n(dir) - dbv.y;
      dK_cu(3,i) = -fe.d2NdX2(k,1,2)*n(dir) - dbv.z;
    }

#if INT_DEBUG > 1
  std::cout <<"dE_cu:"<< dE_cu <<"dK_cu:"<< dK_cu;
#endif

  Matrix dE_ca, dK_ca;
  dE_ca.multiply(T,dE_cu);
  dK_ca.multiply(T,dK_cu);

  // Calculate internal forces
  Vector fiem, fieb;
  dE_ca.multiply(N_ca,fiem,true); // fiem = dE_ca^t * N_ca
  dK_ca.multiply(M_ca,fieb,true); // fieb = dK_ca^t * M_ca
  ES.add(fiem,-fe.detJxW).add(fieb,-fe.detJxW);
#if INT_DEBUG > 1
  std::cout <<"fiem:"<< fiem <<"fieb:"<< fieb;
#endif

  if (EK.empty())
    return true;

  Matrix3D ddE_ca(3,nedof,nedof), ddK_ca(3,nedof,nedof);
  for (size_t kr = 1; kr <= nenod; kr++)
    for (int dirr = 1; dirr <= 3; dirr++)
    {
      size_t r = 3*kr-3 + dirr;
      for (size_t s = 1; s <= r; s++)
      {
        size_t ks = (s-1)/3 + 1;
        int dirs = (s-1)%3 + 1;
        if (dirr == dirs) {
          Vec3 ddE_cu; // Strain
          ddE_cu(1) = fe.dNdX(kr,1)*fe.dNdX(ks,1);
          ddE_cu(2) = fe.dNdX(kr,2)*fe.dNdX(ks,2);
          ddE_cu(3) = 0.5*(fe.dNdX(kr,1)*fe.dNdX(ks,2) +
                           fe.dNdX(kr,2)*fe.dNdX(ks,1));
          ddE_ca.fillColumn(r,s,(T*ddE_cu).vec());
        }

        Vec3 ddg3;
        int dirt = 6-dirr-dirs;
        int ddir = dirr-dirs;
        if (ddir == -1 || ddir == 2)
          ddg3(dirt) =  fe.dNdX(kr,1)*fe.dNdX(ks,2) - fe.dNdX(ks,1)*fe.dNdX(kr,2);
        else if (ddir == 1 || ddir == -2)
          ddg3(dirt) = -fe.dNdX(kr,1)*fe.dNdX(ks,2) + fe.dNdX(ks,1)*fe.dNdX(kr,2);
        double C = -(ddg3*g3 + dg3.getColumn(r)*dg3.getColumn(s))/lg3_3;
        double D = 3.0*g3dg3(r)*g3dg3(s)/lg3_5;
        Vec3 ddn = ddg3/lg3 - (dg3.getColumn(r)*g3dg3lg3_3(s) +
                               g3dg3lg3_3(r)*dg3.getColumn(s)) + (C+D)*g3;
        Vec3 ddbv = ddn * (Hn.empty() ? H0 : Hn);

        Vec3 ddK_cu; // Curvature
        ddK_cu(1) = -fe.d2NdX2(kr,1,1)*dn(dirr,s) - fe.d2NdX2(ks,1,1)*dn(dirs,r) - ddbv.x;
        ddK_cu(2) = -fe.d2NdX2(kr,2,2)*dn(dirr,s) - fe.d2NdX2(ks,2,2)*dn(dirs,r) - ddbv.y;
        ddK_cu(3) = -fe.d2NdX2(kr,1,2)*dn(dirr,s) - fe.d2NdX2(ks,1,2)*dn(dirs,r) - ddbv.z;
        ddK_ca.fillColumn(r,s,(T*ddK_cu).vec());
      }
    }

  // Calculate tangent stiffness matrix
  Matrix dN_ca, dM_ca, kem(nedof,nedof), keb(nedof,nedof);
  dN_ca.multiply(Dm,dE_ca);
  dM_ca.multiply(Db,dK_ca);
  for (size_t r = 1; r <= nedof; r++)
    for (size_t s = 1; s <= r; s++)
    {
      kem(r,s) = N_ca*ddE_ca.getColumn(r,s);
      keb(r,s) = M_ca*ddK_ca.getColumn(r,s);
      if (s < r)
      {
        kem(s,r) = kem(r,s);
        keb(s,r) = keb(r,s);
      }
    }

  kem.multiply(dN_ca,dE_ca,true,false,true); // kem += dN_ca^t * dE_ca
  keb.multiply(dM_ca,dK_ca,true,false,true); // keb += dM_ca^t * dK_ca
  EK.add(kem,fe.detJxW).add(keb,fe.detJxW);

  return true;
}


bool NLKirchhoffLoveShell::evalSol (Vector& sm, Vector& sb, const Vectors& eV,
                                    const FiniteElement& fe, const Vec3& X,
                                    bool) const
{
  if (eV.size() < 2 || eV.front().empty() || eV.back().empty())
  {
    std::cerr <<" *** NLKirchhoffLoveShell::evalSol: No displacement vectors."
              << std::endl;
    return false;
  }
  else if (eV.front().size() != 3*fe.d2NdX2.dim(1))
  {
    std::cerr <<" *** NLKirchhoffLoveShell::evalSol: Invalid solution vector."
              <<"\n     size(eV) = "<< eV.front().size() <<"   size(d2NdX2) = "
              << fe.d2NdX2.dim(1) <<","<< fe.d2NdX2.dim(2)*fe.d2NdX2.dim(3)
              << std::endl;
    return false;
  }

  // Co-variant basis and Hessian in deformed configuration
  Matrix Gd, Hd;
  Matrix3D Hess;
  Gd.multiplyMat(eV.back(),fe.dNdX);
  if (Hess.multiplyMat(eV.back(),fe.d2NdX2))
    utl::Hessian(Hess,Hd);
  else
    return false;

  Matrix Dm, Db;
  if (!this->formDmatrix(Dm,Db,fe,X))
    return false;

  // Calculate the metrics of the initial and deformed configurations
  Vec3 g1, g2, g3, n, gab, Gab; Matrix T(3,3);
  this->getMetrics(fe.G,g1,g2,g3,n,Gab,&T);
  this->getMetrics(Gd,g1,g2,g3,n,gab);
  Vec3 Bv = n * fe.H;
  Vec3 bv = n * Hd;

#if INT_DEBUG > 1
  std::cout <<"\nNLKirchhoffLoveShell::evalSol(X="<< X
            <<")\n\tg1 = "<< g1 <<"\n\tg2 = "<< g2 <<"\n\tg3 = "<< g3
            <<"\n\tn = "<< n <<"\n\tgab = "<< gab <<"\n\tgab = "<< bv <<"\n";
#endif

  // Strain vectors referred to curvilinear coordinate system
  Vec3 E_cu = 0.5*(gab-Gab);
  Vec3 K_cu = (Bv - bv);

  // Strain vectors referred to cartesian coordinate system
  Vec3 E_ca = T * E_cu;
  Vec3 K_ca = T * K_cu;

  // Stress resultants referred to cartesian coordinate system
  if (!Dm.multiply(E_ca.vec(),sm)) // sm = Dm*E_ca
    return false;
  if (!Db.multiply(K_ca.vec(),sb)) // sb = Db*K_ca
    return false;

#if INT_DEBUG > 1
  std::cout <<"\tE_ca = "<< E_ca <<"\n\tK_ca = "<< K_ca
            <<"\nN_ca:"<< sm <<"M_ca:"<< sb;
#endif
  return true;
}
