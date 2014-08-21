// $Id$
//==============================================================================
//!
//! \file LinearElasticity.C
//!
//! \date Nov 12 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for linear elasticity problems.
//!
//==============================================================================

#include "LinearElasticity.h"
#include "MaterialBase.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "Functions.h"
#include "Utilities.h"
#include "VTF.h"
#include "tinyxml.h"


LinearElasticity::LinearElasticity (unsigned short int n, bool axS, bool GPout)
  : Elasticity(n,axS)
{
  myTemp0 = myTemp = NULL;
  myItgPts = n == 2 && GPout ? new Vec3Vec() : NULL;
}


bool LinearElasticity::parse (const TiXmlElement* elem)
{
  bool initT = !strcasecmp(elem->Value(),"initialtemperature");
  if (initT || !strcasecmp(elem->Value(),"temperature"))
  {
    std::string type;
    utl::getAttribute(elem,"type",type,true);
    const TiXmlNode* tval = elem->FirstChild();
    if (tval)
    {
      if (initT)
      {
        std::cout <<"\tInitial temperature";
        myTemp0 = utl::parseRealFunc(tval->Value(),type);
      }
      else
      {
        std::cout <<"\tTemperature";
        myTemp = utl::parseRealFunc(tval->Value(),type);
      }
      std::cout << std::endl;
    }
  }

  return true;
}


void LinearElasticity::setMode (SIM::SolutionMode mode)
{
  if (mode == SIM::RECOVERY && m_mode != mode)
  {
    maxVal.resize(this->getNoFields(2));
    std::fill(maxVal.begin(),maxVal.end(),PointValue(Vec3(),0.0));
  }

  this->ElasticBase::setMode(mode);

  // These quantities are not needed in linear problems
  if (mode != SIM::BUCKLING) eKg = 0;
  if (mode != SIM::DYNAMIC)  iS  = 0;
}


void LinearElasticity::initIntegration (size_t nGp, size_t nBp)
{
  this->Elasticity::initIntegration(nGp,nBp);
  if (myItgPts) myItgPts->resize(nGp);
}


bool LinearElasticity::hasTractionValues() const
{
  if (myItgPts && !myItgPts->empty())
    return true;

  return this->Elasticity::hasTractionValues();
}


bool LinearElasticity::writeGlvT (VTF* vtf, int iStep,
                                  int& geoBlk, int& nBlock) const
{
  bool ok = this->Elasticity::writeGlvT(vtf,iStep,geoBlk,nBlock);
  if (ok && vtf && myItgPts && !myItgPts->empty())
    return vtf->writePoints(*myItgPts,geoBlk);

  return ok;
}


bool LinearElasticity::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                                const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  bool lHaveStrains = false;
  SymmTensor eps(nsd,axiSymmetry), sigma(nsd,axiSymmetry);

  Matrix Bmat, Cmat;
  if (eKm || eKg || iS || (eS && myTemp))
  {
    // Compute the strain-displacement matrix B from N, dNdX and r = X.x,
    // and evaluate the symmetric strain tensor if displacements are available
    if (!this->kinematics(elMat.vec.front(),fe.N,fe.dNdX,X.x,Bmat,eps,eps))
      return false;
    else if (!eps.isZero(1.0e-16))
      lHaveStrains = true;

    // Evaluate the constitutive matrix and the stress tensor at this point
    double U;
    if (!material->evaluate(Cmat,sigma,U,fe,X,eps,eps))
      return false;

#if INT_DEBUG > 3
    std::cout <<"LinearElasticity::evalInt(X = "<< X <<")\n"
              <<"Bmat ="<< Bmat <<"Cmat ="<< Cmat;
#if INT_DEBUG > 4
    if (lHaveStrains) std::cout <<"sigma =\n"<< sigma;
#endif
    std::cout << std::endl;
#endif
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
  {
    // Integrate the geometric stiffness matrix
    double r = axiSymmetry ? X.x + elMat.vec.front().dot(fe.N,0,nsd) : 0.0;
    this->formKG(elMat.A[eKg-1],fe.N,fe.dNdX,r,sigma,detJW);
  }

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
  {
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(elMat.b[eS-1],fe.N,X,detJW);
    // Integrate the load vector due to initial or temperature strains
    if (!this->formInitStrainForces(elMat,fe.N,Bmat,Cmat,X,detJW))
      return false;
  }

  // Store Gauss point coordinates for visualization
  if (myItgPts && fe.iGP < myItgPts->size())
    myItgPts->at(fe.iGP) = X;

  return true;
}


bool LinearElasticity::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                                const Vec3& X, const Vec3& normal) const
{
  if (!tracFld && !fluxFld)
  {
    std::cerr <<" *** LinearElasticity::evalBou: No tractions."<< std::endl;
    return false;
  }
  else if (!eS)
  {
    std::cerr <<" *** LinearElasticity::evalBou: No load vector."<< std::endl;
    return false;
  }

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  const double detJW = axiSymmetry ? 2.0*M_PI*X.x*fe.detJxW : fe.detJxW;

  // Evaluate the surface traction
  Vec3 T = this->getTraction(X,normal);

  // Store traction value for visualization
  if (fe.iGP < tracVal.size() && !T.isZero())
  {
    tracVal[fe.iGP].first = X;
    tracVal[fe.iGP].second += T;
  }

  // Integrate the force vector
  Vector& ES = static_cast<ElmMats&>(elmInt).b[eS-1];
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (unsigned short int i = 1; i <= nsd; i++)
      ES(nsd*(a-1)+i) += T[i-1]*fe.N(a)*detJW;

  return true;
}


/*!
  This method evaluates the stabilization term used in immersed boundary
  simulations. According to Mats Larsons suggestion.
*/

bool LinearElasticity::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                                const Vec3& X, const Vec3&) const
{
  if (!eKm)
  {
    std::cerr <<" *** LinearElasticity::evalInt: No material stiffness matrix."
              << std::endl;
    return false;
  }
  else if (axiSymmetry)
  {
    std::cerr <<" *** LinearElasticity::evalInt: Axi-symmetric problem!"
              << std::endl;
    return false;
  }

  int pdir = 0;
  if (fe.xi == -1.0)
    pdir = -1;
  else if (fe.xi == 1.0)
    pdir =  1;
  else if (fe.eta == -1.0)
    pdir = -2;
  else if (fe.eta ==  1.0)
    pdir =  2;
  else if (fe.zeta == -1.0)
    pdir = -3;
  else if (fe.zeta ==  1.0)
    pdir =  3;
  else
  {
    std::cerr <<" *** LinearElasticity::evalInt: Not on an interface, "
              <<" xi="<< fe.xi <<" eta="<< fe.eta <<" zeta="<< fe.zeta
              << std::endl;
    return false;
  }

  // Compute the element length in the parametric normal direction
  // of the interface, assuming here a Cartesian grid, for now...
  double h = 0.0;
  if (pdir == 1 || pdir == -1)
    h = (fe.XC[1] - fe.XC.front()).length();
  else if (pdir == 2 || pdir == -2)
    h = (fe.XC[2] - fe.XC.front()).length();
  else
    h = (fe.XC[4] - fe.XC.front()).length();

  // Evaluate the stiffness at current point
  double E = material->getStiffness(X)*gamma;

  // Evaluate the stabilization constant
  double hJW = pow(h,2*fe.p+1)*E*fe.detJxW;
  if (pdir < 0) hJW = -hJW;

  // Integrate the interface jump term
  Matrix& EK = static_cast<ElmMats&>(elmInt).A[eKm-1];
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (size_t b = 1; b <= fe.N.size(); b++)
      for (unsigned short int i = 1; i <= nsd; i++)
        EK(nsd*(a-1)+i,nsd*(b-1)+i) += hJW*fe.N(a)*fe.N(b);

  return true;
}


int LinearElasticity::getIntegrandType () const
{
  return INTERFACE_TERMS | ELEMENT_CORNERS | NORMAL_DERIVS;
}


double LinearElasticity::getThermalStrain (const Vector&, const Vector&,
                                           const Vec3& X) const
{
  if (!myTemp) return 0.0;

  double T0 = myTemp0 ? (*myTemp0)(X) : 0.0;
  double T = (*myTemp)(X);
  return material->getThermalExpansion(T)*(T-T0);
}


bool LinearElasticity::formInitStrainForces (ElmMats& elMat, const Vector& N,
                                             const Matrix& B, const Matrix& C,
                                             const Vec3& X, double detJW) const
{
  if (!eS || !myTemp)
    return true; // No temperature field

  // Strains due to thermal expansion
  SymmTensor eps(nsd,axiSymmetry);
  eps = this->getThermalStrain(N,N,X)*detJW;

  // Stresses due to thermal expansion
  Vector sigma0;
  if (!C.multiply(eps,sigma0))
    return false;

  // Integrate external forces due to thermal expansion
  SymmTensor sigma(nsd,axiSymmetry); sigma = sigma0;
  return B.multiply(sigma,elMat.b[eS-1],true,true); // ES += B^T*sigma0
}
