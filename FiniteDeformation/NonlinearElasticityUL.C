// $Id$
//==============================================================================
//!
//! \file NonlinearElasticityUL.C
//!
//! \date Sep 21 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity problems.
//!
//==============================================================================

#include "NonlinearElasticityUL.h"
#include "MaterialBase.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "TimeDomain.h"
#include "Tensor.h"
#include "Function.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include "int_debug.h"

#ifndef epsR
//! \brief Zero tolerance for the radial coordinate.
#define epsR 1.0e-16
#endif


NonlinearElasticityUL::NonlinearElasticityUL (unsigned short int n,
					      bool axS, char lop)
  : Elasticity(n,axS), loadOp(lop)
{
  // Only the current solution is needed
  primsol.resize(1);
}


void NonlinearElasticityUL::printLog () const
{
  IFEM::cout <<"NonlinearElasticityUL: Updated Lagrangian formulation"
             << std::endl;

  this->Elasticity::printLog();
}


void NonlinearElasticityUL::setMode (SIM::SolutionMode mode)
{
  if (mode == SIM::RECOVERY && mode != m_mode)
    this->initMaxVals();

  if (mode == SIM::RHS_ONLY)
  {
    this->Elasticity::setMode(SIM::STATIC);
    m_mode = SIM::RHS_ONLY;
    primsol.resize(1);
  }
  else
    this->Elasticity::setMode(mode);
}


void NonlinearElasticityUL::initIntegration (size_t nGp, size_t nBp)
{
  if (material)
    material->initIntegration(nGp);

  this->Elasticity::initIntegration(nGp,nBp);
}


void NonlinearElasticityUL::initIntegration (const TimeDomain& prm,
					     const Vector&, bool)
{
  if (material)
    material->initIntegration(prm);
}


void NonlinearElasticityUL::initResultPoints (double lambda, char prinDir)
{
  if (material && prinDir >= 0)
    material->initResultPoints();

  this->Elasticity::initResultPoints(lambda,prinDir);
}


bool NonlinearElasticityUL::evalInt (LocalIntegral& elmInt,
				     const FiniteElement& fe,
				     const TimeDomain& prm,
				     const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  // Evaluate the deformation gradient, F, and the Green-Lagrange strains, E
  Matrix Bmat, dNdx;
  Tensor F(nDF);
  SymmTensor E(nsd,axiSymmetry);
  if (!this->kinematics(elMat.vec.front(),fe.N,fe.dNdX,X.x,Bmat,F,E))
    return false;

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  double detJW = axiSymmetry ? 2.0*M_PI*X.x*fe.detJxW : fe.detJxW;
  double r = axiSymmetry ? X.x : 0.0;

  bool lHaveStrains = !E.isZero(1.0e-16);
  if (lHaveStrains)
  {
    // Invert the deformation gradient ==> Fi
    Matrix Fi(nsd,nsd);
    if (nDF == nsd)
      Fi.fill(F.ptr());
    else
      for (unsigned short int i = 1; i <= nsd; i++)
        for (unsigned short int j = 1; j <= nsd; j++)
          Fi(i,j) = F(i,j);

    double J = Fi.inverse();
    if (axiSymmetry) J *= F(3,3);
    if (J == 0.0) return false;

    // Scale with J=|F| since we are integrating over current configuration
    detJW *= J;

    if (eKm || iS)
    {
      // Push-forward the basis function gradients to current configuration
      dNdx.multiply(fe.dNdX,Fi); // dNdx = dNdX * F^-1
      // Compute the small-deformation strain-displacement matrix B from dNdx
      if (axiSymmetry)
      {
	r += elMat.vec.front().dot(fe.N,0,nsd);
	this->formBmatrix(Bmat,fe.N,dNdx,r);
      }
      else
	this->formBmatrix(Bmat,dNdx);

#if INT_DEBUG > 0
      std::cout <<"NonlinearElasticityUL::dNdx ="<< dNdx;
      std::cout <<"NonlinearElasticityUL::B ="<< Bmat;
#endif
    }
  }
  else if (eKm || iS)
  {
    // Initial state, no deformation yet
    if (axiSymmetry)
      this->formBmatrix(Bmat,fe.N,fe.dNdX,r);
    else
      this->formBmatrix(Bmat,fe.dNdX);
  }

  // Evaluate the constitutive relation
  Matrix Cmat;
  SymmTensor sigma(nsd,axiSymmetry);
  if (eKm || eKg || iS)
  {
    double U = 0.0;
    if (!material->evaluate(Cmat,sigma,U,fe,X,F,E,(eKg || iS),&prm))
      return false;
  }

  if (eKm)
  {
    // Integrate the material stiffness matrix
    Matrix CB;
    CB.multiply(Cmat,Bmat).multiply(detJW); // CB = C*B*|J|*w
    elMat.A[eKm-1].multiply(Bmat,CB,true,false,true); // EK += B^T * CB
  }

  if (eKg && lHaveStrains)
    // Integrate the geometric stiffness matrix
    this->formKG(elMat.A[eKg-1],fe.N,dNdx,r,sigma,detJW);

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(elMat.A[eM-1],fe.N,X,detJW);

  if (iS && lHaveStrains)
  {
    // Integrate the internal forces
    sigma *= -detJW;
    if (!Bmat.multiply(sigma,elMat.b[iS-1],true,true)) // ES -= B^T*sigma
      return false;
  }

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(elMat.b[eS-1],elMat.c,fe.N,X,detJW);

  if (gS)
    // Integrate the load gradient vector due to other body forces
    this->formBodyForce(elMat.b[eS-1],elMat.c,fe.N,X,detJW,true);

  return true;
}


bool NonlinearElasticityUL::evalBou (LocalIntegral& elmInt,
				     const FiniteElement& fe,
				     const Vec3& X, const Vec3& normal) const
{
  if (!tracFld && !fluxFld)
  {
    std::cerr <<" *** NonlinearElasticityUL::evalBou: No tractions."
	      << std::endl;
    return false;
  }
  else if (!eS)
  {
    std::cerr <<" *** NonlinearElasticityUL::evalBou: No load vector."
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

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  double detJW = axiSymmetry ? 2.0*M_PI*X.x*fe.detJxW : fe.detJxW;

  if (loadOp == 1)
  {
    // Compute the deformation gradient, F
    Tensor F(nDF);
    if (!this->formDefGradient(elmInt.vec.front(),fe.N,fe.dNdX,X.x,F))
      return false;

    // Check for with-rotated pressure load
    if (tracFld && tracFld->isNormalPressure())
    {
      // Compute its inverse and determinant, J
      double J = F.inverse();
      if (J == 0.0) return false;

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
    else
      // Scale with J=|F| since we are integrating over current configuration
      detJW *= F.det();
  }

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


bool NonlinearElasticityUL::kinematics (const Vector& eV,
					const Vector& N, const Matrix& dNdX,
					double r, Matrix&, Tensor& F,
					SymmTensor& E) const
{
  // Compute the deformation gradient, [F] = [I] + [dudX] = [I] + [dNdX]*[u],
  if (eV.empty())
  {
    // Initial state, unit deformation gradient and zero strains
    F = 1.0;
    E.zero();
    return true;
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
  std::cout <<"NonlinearElasticityUL::F =\n"<< F;
#endif

  return true;
}


bool NonlinearElasticityUL::diverged (size_t iP) const
{
  return material ? material->diverged(iP) : false;
}


NormBase* NonlinearElasticityUL::getNormIntegrand (AnaSol*) const
{
  return new ElasticityNormUL(const_cast<NonlinearElasticityUL&>(*this));
}


void ElasticityNormUL::initIntegration (size_t nGp, size_t nBp)
{
  Ux.resize(nBp,0.0);
  up.resize(nBp);
  tp.resize(nBp);

  this->ElasticityNorm::initIntegration(nGp,nBp);
}


bool ElasticityNormUL::evalInt (LocalIntegral& elmInt,
				const FiniteElement& fe,
				const TimeDomain& prm,
				const Vec3& X) const
{
  NonlinearElasticityUL& ulp = static_cast<NonlinearElasticityUL&>(myProblem);

  // Evaluate the deformation gradient, F, and the Green-Lagrange strains, E
  Matrix B;
  Tensor F(ulp.nDF);
  SymmTensor E(ulp.nDF);
  if (!ulp.kinematics(elmInt.vec.front(),fe.N,fe.dNdX,X.x,B,F,E))
    return false;

  // Compute the strain energy density, U(E) = Int_E (S:Eps) dEps
  // and the Cauchy stress tensor, sigma
  Matrix Cmat; double U = 0.0;
  SymmTensor sigma(E.dim(),ulp.isAxiSymmetric()||ulp.material->isPlaneStrain());
  if (!ulp.material->evaluate(Cmat,sigma,U,fe,X,F,E,3,&prm))
    if (!ulp.material->diverged(fe.iGP+1))
      return false;

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  double detJW = ulp.isAxiSymmetric() ? 2.0*M_PI*X.x*fe.detJxW : fe.detJxW;

  // Integrate the norms
  return evalInt(static_cast<ElmNorm&>(elmInt),sigma,U,F.det(),detJW);
}


size_t ElasticityNormUL::getNoFields (int group) const
{
  return group < 1 ? 1 : 7;
}


bool ElasticityNormUL::evalInt (ElmNorm& pnorm, const SymmTensor& S,
				double U, double detF, double detJxW)
{
  // Integrate the energy norm a(u^h,u^h) = Int_Omega0 U(E) dV0
  pnorm[0] += U*detJxW;

  // Integrate the L2-norm ||S|| = Int_Omega S:S dV
  detJxW *= detF;
  pnorm[2] += S.L2norm(false)*detJxW;

  // Integrate the L2-norm ||p|| = Int_Omega (trace(S)/nsd)^2 dV
  double p = S.trace() / (double)(S.size() == 4 ? 3 : S.dim());
  pnorm[3] += p*p*detJxW;

  // Integrate the L2-norm ||S_dev|| = Int_Omega (S-p*I):(S-p*I) dV
  SymmTensor Sdev(S); Sdev -= p;
  pnorm[4] += Sdev.L2norm(false)*detJxW;

  // Integrate the von Mises stress norm
  pnorm[5] += S.vonMises(false)*detJxW;

  // Integrate the volume
  pnorm[6] += detJxW;

  return true;
}


std::string ElasticityNormUL::getName (size_t i, size_t j,
                                       const char* prefix) const
{
  if (i != 1 || j == 0 || j > 7)
    return this->ElasticityNorm::getName(i,j,prefix);

  static const char* s[7] = {
    "a(u^h,u^h)^0.5",
    "((f,u^h)+(t,u^h))^0.5",
    "(s^h,s^h)^0.5",
    "(p^h,p^h)^0.5",
    "(e^h,e^h), e^h=S^h-p^h*I",
    "vm(s^h)",
    "volume"
  };

  if (!prefix)
    return s[j-1];

  return prefix + std::string(" ") + s[j-1];
}


bool ElasticityNormUL::evalBou (LocalIntegral& elmInt,
				const FiniteElement& fe,
				const Vec3& X, const Vec3& normal) const
{
  NonlinearElasticityUL& ulp = static_cast<NonlinearElasticityUL&>(myProblem);
  if (!ulp.haveLoads()) return true;

  // Evaluate the current surface traction
  Vec3 t = ulp.getTraction(X,normal);
  // Evaluate the current displacement field
  Vec3 u = ulp.evalSol(elmInt.vec.front(),fe.N);

  // Integrate the external energy (path integral)
  size_t iP = fe.iGP;
#ifdef INDEX_CHECK
  if (iP >= Ux.size())
  {
    std::cerr <<" *** ElasticityNormUL::evalBou: Integration point "<< iP+1
	      <<" out of range [1,"<< Ux.size() <<"]."<< std::endl;
    return false;
  }
#endif
  Ux[iP] += 0.5*(t+tp[iP])*(u-up[iP]);
  tp[iP] = t;
  up[iP] = u;

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  double detJW = ulp.isAxiSymmetric() ? 2.0*M_PI*X.x*fe.detJxW : fe.detJxW;

  static_cast<ElmNorm&>(elmInt)[1] += Ux[iP]*detJW;
  return true;
}
