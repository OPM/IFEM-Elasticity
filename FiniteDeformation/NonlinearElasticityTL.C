// $Id$
//==============================================================================
//!
//! \file NonlinearElasticityTL.C
//!
//! \date May 25 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity problems.
//!
//==============================================================================

#include "NonlinearElasticityTL.h"
#include "MaterialBase.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "Tensor.h"
#include "Function.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include "int_debug.h"

#ifndef epsR
//! \brief Zero tolerance for the radial coordinate.
#define epsR 1.0e-16
#endif


NonlinearElasticityTL::NonlinearElasticityTL (unsigned short int n, bool axS,
                                              bool evalAtElmCenter)
  : Elasticity(n,axS)
{
  myIntegrandType = evalAtElmCenter ? ELEMENT_CENTER : STANDARD;
  formB = false;
}


void NonlinearElasticityTL::printLog () const
{
  IFEM::cout <<"NonlinearElasticityTL: Total Lagrangian formulation"
             << std::endl;

  this->Elasticity::printLog();
}


void NonlinearElasticityTL::setMode (SIM::SolutionMode mode)
{
  if (mode == SIM::RHS_ONLY)
  {
    this->Elasticity::setMode(SIM::STATIC);
    m_mode = SIM::RHS_ONLY;
    primsol.resize(1);
  }
  else
    this->Elasticity::setMode(mode);

  switch (mode)
    {
    case SIM::ARCLEN:
    case SIM::STATIC:
    case SIM::DYNAMIC:
    case SIM::BUCKLING:
    case SIM::RHS_ONLY:
      formB = true;
      break;

    default:
      formB = false;
    }
}


bool NonlinearElasticityTL::evalInt (LocalIntegral& elmInt,
				     const FiniteElement& fe,
				     const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  // Evaluate the deformation gradient, F, and the Green-Lagrange strains, E,
  // and compute the nonlinear strain-displacement matrix, B, from dNdX and F
  Matrix Bmat;
  Tensor F(nDF);
  SymmTensor E(nsd,axiSymmetry), S(nsd,axiSymmetry);
  if (!this->kinematics(elMat.vec.front(),fe.N,fe.dNdX,X.x,Bmat,F,E))
    return false;

  // Evaluate the constitutive relation
  Matrix Cmat;
  bool lHaveStrains = !E.isZero(1.0e-16);
  if (eKm || eKg || iS)
  {
    double U;
    int iopm = lHaveStrains && (eKg || iS) ? 2 : 0;
    if (myIntegrandType & ELEMENT_CENTER && fe.XC.size() == 2)
    {
      // Evaluate material properties at the element center.
      // The Cartesian coordinates are assumed store in fe.XC(1), and
      // the parametric coordinates w.r.t. current spline patch in fe.XC(2)
      FiniteElement fe0(0,fe.iGP);
      fe0.u = fe.XC.back().x;
      fe0.v = fe.XC.back().y;
      fe0.w = fe.XC.back().z;
      Vec4 X0(fe.XC.front());
      X0.t = static_cast<const Vec4&>(X).t;
      if (!material->evaluate(Cmat,S,U,fe0,X0,F,E,iopm))
        return false;
    }
    else if (!material->evaluate(Cmat,S,U,fe,X,F,E,iopm))
      return false;
  }

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  const double detJW = axiSymmetry ? 2.0*M_PI*X.x*fe.detJxW : fe.detJxW;

  if (eKm)
  {
    // Integrate the material stiffness matrix
    Matrix CB;
    CB.multiply(Cmat,Bmat).multiply(detJW); // CB = C*B*|J|*w
    elMat.A[eKm-1].multiply(Bmat,CB,true,false,true); // EK += B^T * CB
  }

  if (eKg && lHaveStrains)
    // Integrate the geometric stiffness matrix
    this->formKG(elMat.A[eKg-1],fe.N,fe.dNdX,X.x,S,detJW);

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(elMat.A[eM-1],fe.N,X,detJW);

  if (iS && lHaveStrains)
  {
    // Integrate the internal forces
    S *= -detJW;
    if (!Bmat.multiply(S,elMat.b[iS-1],true,true)) // ES -= B^T*S
      return false;
  }

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(elMat.b[eS-1],elMat.c,fe.N,X,detJW);

  if (gS)
    // Integrate the load gradient vector due to other body forces
    this->formBodyForce(elMat.b[gS-1],elMat.c,fe.N,X,detJW,true);

  return true;
}


bool NonlinearElasticityTL::evalBou (LocalIntegral& elmInt,
				     const FiniteElement& fe,
				     const Vec3& X, const Vec3& normal) const
{
  if (!tracFld && !fluxFld)
  {
    std::cerr <<" *** NonlinearElasticityTL::evalBou: No tractions."
	      << std::endl;
    return false;
  }
  else if (!eS)
  {
    std::cerr <<" *** NonlinearElasticityTL::evalBou: No load vector."
	      << std::endl;
    return false;
  }

  // Evaluate the surface traction
  Vec3 Tg, T = this->getTraction(X,normal);
  if (gS) Tg = this->getTraction(X,normal,true);

  // Store traction value for visualization
  if (fe.iGP < tracVal.size() && !T.isZero())
  {
    tracVal[fe.iGP].first = X;
    tracVal[fe.iGP].second += T;
  }

  // Check for with-rotated pressure load
  if (tracFld && tracFld->isNormalPressure())
  {
    // Compute the deformation gradient, F
    Tensor F(nDF);
    if (!this->formDefGradient(elmInt.vec.front(),fe.N,fe.dNdX,X.x,F))
      return false;

    // Compute its inverse and determinant, J
    double J = F.inverse();
    if (J <= 0.0) return false;

    // Pull-back the normal traction to the initial configuration.
    // See equation (3.4.5) on page 102 in T. Belytschko et. al (2000):
    //    p*n*dS ==> J*p*n0*F^-1*dS0
    auto&& pullBack=[F,J,this](Vec3& T)
    {
      Vec3 t = J*T; T = 0.0;
      for (unsigned short int i = 1; i <= nsd; i++)
        for (unsigned short int j = 1; j <= nsd; j++)
          T[i-1] += t[j-1]*F(j,i);
    };
    pullBack(T);
    if (gS) pullBack(Tg);
  }

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  const double detJW = axiSymmetry ? 2.0*M_PI*X.x*fe.detJxW : fe.detJxW;

  // Integrate the force vector
  Vector& ES = static_cast<ElmMats&>(elmInt).b[eS-1];
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (unsigned short int i = 1; i <= nsd; i++)
      ES(nsd*(a-1)+i) += T[i-1]*fe.N(a)*detJW;

  // Integrate total external load
  RealArray& sumLoad = static_cast<ElmMats&>(elmInt).c;
  for (unsigned short int i = 0; i < nsd && i < sumLoad.size(); i++)
    sumLoad[i] += T[i]*detJW;

  if (gS)
  {
    // Integrate the gradient force vector
    Vector& GS = static_cast<ElmMats&>(elmInt).b[gS-1];
    for (size_t a = 1; a <= fe.N.size(); a++)
      for (unsigned short int i = 1; i <= nsd; i++)
        GS(nsd*(a-1)+i) += Tg[i-1]*fe.N(a)*detJW;
  }

  return true;
}


bool NonlinearElasticityTL::kinematics (const Vector& eV,
					const Vector& N, const Matrix& dNdX,
					double r, Matrix& Bmat, Tensor& F,
					SymmTensor& E) const
{
  // Compute the deformation gradient, [F] = [I] + [dudX] = [I] + [dNdX]*[u],
  if (eV.empty())
  {
    // Initial state, unit deformation gradient and linear B-matrix
    F = 1.0;
    E.zero();
    if (axiSymmetry)
      return this->formBmatrix(Bmat,N,dNdX,r);
    else
      return this->formBmatrix(Bmat,dNdX);
  }
  else if (!this->formDefGradient(eV,N,dNdX,r,F,true))
    return false;

  // Form the Green-Lagrange strain tensor, E_ij = 0.5*(F_ij+F_ji+F_ki*F_kj).
  // Note that for the shear terms (i/=j) we actually compute 2*E_ij
  // to be consistent with the engineering strain style constitutive matrix.
  // TODO: How is this for axisymmetric problems?
  unsigned short int i, j, k;
  for (i = 1; i <= E.dim(); i++)
    for (j = 1; j <= i; j++)
    {
      double Eij = F(i,j) + F(j,i);
      for (k = 1; k <= nsd; k++)
        Eij += F(k,i)*F(k,j);
      E(i,j) = i == j ? 0.5*Eij : Eij;
    }

  // Add the unit tensor to F to form the deformation gradient
  F += 1.0;
  // Add the dU/r term to the F(3,3)-term for axisymmetric problems
  if (axiSymmetry && r > epsR) F(3,3) += eV.dot(N,0,nsd)/r;

#if INT_DEBUG > 0
  std::cout <<"NonlinearElasticityTL::F =\n"<< F;
#endif

  if (!formB || E.dim() < nsd) return true;

  // Form the nonlinear B-matrix
  const size_t nenod = dNdX.rows();
  const size_t nstrc = axiSymmetry ? 4 : nsd*(nsd+1)/2;
  Bmat.resize(nstrc*nsd,nenod,true);

#define INDEX(i,j) i+nstrc*(j-1)

  // Normal strain part
  for (size_t a = 1; a <= nenod; a++)
    for (i = 1; i <= nsd; i++)
      for (j = 1; j <= nsd; j++)
	Bmat(INDEX(j,i),a) = F(i,j)*dNdX(a,j);

  // Shear strain part
  if (nsd == 3)
    for (size_t a = 1; a <= nenod; a++)
      for (i = 1; i <= nsd; i++)
      {
	Bmat(INDEX(4,i),a) = F(i,1)*dNdX(a,2) + F(i,2)*dNdX(a,1);
	Bmat(INDEX(5,i),a) = F(i,2)*dNdX(a,3) + F(i,3)*dNdX(a,2);
	Bmat(INDEX(6,i),a) = F(i,3)*dNdX(a,1) + F(i,1)*dNdX(a,3);
      }

  else if (nsd == 2)
  {
    for (size_t a = 1; a <= nenod; a++)
      if (axiSymmetry)
      {
	for (i = 1; i <= nsd; i++)
	  Bmat(INDEX(4,i),a) = F(i,1)*dNdX(a,2) + F(i,2)*dNdX(a,1);
	// Hoop strain part for axisymmetry (TODO: check this)
	Bmat(INDEX(3,1),i) = F(3,3) * (r <= epsR ? dNdX(a,1) : N(a)/r);
      }
      else
	for (i = 1; i <= nsd; i++)
	  Bmat(INDEX(3,i),a) = F(i,1)*dNdX(a,2) + F(i,2)*dNdX(a,1);
  }

  Bmat.resize(nstrc,nsd*nenod);
#if INT_DEBUG > 0
  std::cout <<"NonlinearElasticityTL::B ="<< Bmat;
#endif
  return true;
}
