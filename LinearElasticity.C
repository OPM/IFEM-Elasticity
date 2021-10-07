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
#include "IFEM.h"
#include "tinyxml.h"


LinearElasticity::LinearElasticity (unsigned short int n, bool axSym,
                                    bool GPout, bool modal)
  : Elasticity(n,axSym)
{
  myTemp0  = myTemp = nullptr;
  myReacIt = nullptr;
  myItgPts = n == 2 && GPout ? new Vec3Vec() : nullptr;
  isModal  = modal;
}


bool LinearElasticity::parse (const TiXmlElement* elem)
{
  bool initT = !strcasecmp(elem->Value(),"initialtemperature");
  if (!initT && strcasecmp(elem->Value(),"temperature"))
    return this->Elasticity::parse(elem);

  std::string type;
  utl::getAttribute(elem,"type",type,true);
  const TiXmlNode* tval = elem->FirstChild();
  if (!tval) return true;

  if (initT)
  {
    IFEM::cout <<"\tInitial temperature";
    myTemp0 = utl::parseRealFunc(tval->Value(),type);
  }
  else
  {
    IFEM::cout <<"\tTemperature";
    myTemp = utl::parseRealFunc(tval->Value(),type);
  }
  IFEM::cout << std::endl;

  return true;
}


void LinearElasticity::setMode (SIM::SolutionMode mode)
{
  if (isModal && mode == SIM::DYNAMIC)
    mode = SIM::RHS_ONLY;

  if (mode >= SIM::RECOVERY && m_mode != mode)
    this->initMaxVals();

  this->ElasticBase::setMode(mode);

  if (dS > 0) dS = mode == SIM::STATIC || mode == SIM::NORMS ? 2 : 1;

  // These quantities are not needed in linear problems
  if (mode != SIM::BUCKLING) eKg = 0;
  if (mode == SIM::STATIC || (isModal && mode == SIM::RHS_ONLY)) iS  = 0;
}


void LinearElasticity::initLHSbuffers (size_t nEl)
{
  if (nEl > 1)
  {
    myKmats.resize(nEl);
    myMmats.resize(nEl);
  }
  else if (nEl == 0 && !myKmats.empty())
  {
    if (eKm > 0) eKm = -eKm;
    if (eKg > 0) eKg = -eKg;
    if (eM  > 0) eM  = -eM;
  }
}


void LinearElasticity::initIntegration (size_t nGp, size_t nBp)
{
  this->Elasticity::initIntegration(nGp,nBp);
  if (myItgPts) myItgPts->resize(nGp);
}


bool LinearElasticity::initElement (const std::vector<int>& MNPC,
                                    const FiniteElement& fe, const Vec3& XC,
                                    size_t, LocalIntegral& elmInt)
{
  if (fe.iel > 0)
  {
    size_t iel = fe.iel - 1;
    if (iel < myKmats.size() && eKm < 0)
      static_cast<ElmMats&>(elmInt).A[-eKm-1] = myKmats[iel];
    if (iel < myMmats.size() && eKm < 0)
      static_cast<ElmMats&>(elmInt).A[-eM-1]  = myMmats[iel];
  }

  size_t nsol = primsol.size();
  while (nsol > 1 && primsol[nsol-1].empty()) nsol--;
  if (nsol <= 1)
    return this->initElement1(MNPC,elmInt.vec);

  int ierr = 0;
  elmInt.vec.resize(nsol);
  for (size_t i = 0; i < nsol && ierr == 0; i++)
    if (!primsol[i].empty())
    {
      bool haveValues = true;
      if (i == 1 && m_mode == SIM::STATIC)
        haveValues = (dualRHS && dualRHS->inDomain(XC));
      else if (i > 0 && m_mode >= SIM::RECOVERY)
        haveValues = (i <= dualFld.size() && dualFld[i-1]->inDomain(XC));
      if (haveValues)
        ierr = utl::gather(MNPC,npv,primsol[i],elmInt.vec[i]);
#if INT_DEBUG > 2
      std::cout <<"Element solution vector "<< i+1 << elmInt.vec[i];
#endif
    }

  if (ierr == 0) return true;

  std::cerr <<" *** LinearElasticity::initElement: Detected "
            << ierr <<" node numbers out of range."<< std::endl;
  return false;
}


GlobalIntegral& LinearElasticity::getGlobalInt (GlobalIntegral* gq) const
{
  if (m_mode == SIM::RHS_ONLY && myReacIt)
    return *myReacIt;

  return this->Elasticity::getGlobalInt(gq);
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
    ok = vtf->writePoints(*myItgPts,geoBlk);

  return ok;
}


bool LinearElasticity::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                                const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  const Vector& eV = elMat.vec.front();

  bool lHaveStrains = false;
  SymmTensor eps(nsd,axiSymmetry), sigma(nsd,axiSymmetry);

  double U = 0.0;
  Matrix Bmat, Cmat;
  if (eKm > 0 || eKg > 0 || (iS > 0 && !eV.empty()) || (eS > 0 && myTemp))
  {
    // Compute the strain-displacement matrix B from N, dNdX and r = X.x,
    // and evaluate the symmetric strain tensor if displacements are available
    if (!this->kinematics(eV,fe.N,fe.dNdX,X.x,Bmat,eps,eps))
      return false;
    else if (!eps.isZero(1.0e-16))
      lHaveStrains = true;

    // Evaluate the constitutive matrix and the stress tensor at this point
    if (!material->evaluate(Cmat,sigma,U,fe,X,eps,eps))
      return false;

#if INT_DEBUG > 3
    std::cout <<"LinearElasticity::evalInt(X = "<< X <<")\nBmat ="<< Bmat;
    if (lHaveStrains) std::cout <<"eps =\n"<< eps;
#if INT_DEBUG > 4
    if (lHaveStrains) std::cout <<"sigma =\n"<< sigma;
#endif
    std::cout <<"Cmat ="<< Cmat << std::endl;
#endif
  }

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  const double detJW = axiSymmetry ? 2.0*M_PI*X.x*fe.detJxW : fe.detJxW;

  if (eKm > 0)
  {
    // Integrate the material stiffness matrix
    Matrix CB;
    CB.multiply(Cmat,Bmat).multiply(detJW); // CB = C*B*|J|*w
    elMat.A[eKm-1].multiply(Bmat,CB,true,false,true); // EK += B^T * CB
  }

  if (eKg > 0 && lHaveStrains)
  {
    // Integrate the geometric stiffness matrix
    double r = axiSymmetry ? X.x + eV.dot(fe.N,0,nsd) : 0.0;
    this->formKG(elMat.A[eKg-1],fe.N,fe.dNdX,r,sigma,detJW);
  }

  if (eM > 0)
    // Integrate the mass matrix
    this->formMassMatrix(elMat.A[eM-1],fe.N,X,detJW);

  if (iS > 0 && lHaveStrains)
  {
    // Integrate the internal forces
    sigma *= -detJW;
    if (!Bmat.multiply(sigma,elMat.b[iS-1],true,true)) // ES -= B^T*sigma
      return false;
  }

  if (eS > 0)
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

  if (dS <= 1 || elMat.vec.size() <= 1)
    return true;

  // Calculate the dual strains and stresses
  if (!this->kinematics(elMat.vec[1],fe.N,fe.dNdX,X.x,Bmat,eps,eps))
    return false;
  else if (eps.isZero(1.0e-16))
    return true; // the extraction function is zero in this element
  else if (!material->evaluate(Cmat,sigma,U,fe,X,eps,eps))
    return false;

#if INT_DEBUG > 3
  std::cout <<"Dual eps =\n"<< eps <<"Dual sigma =\n"<< sigma;
#endif

  // Integrate the dual load vector
  sigma *= detJW;
  return Bmat.multiply(sigma,elMat.b[dS-1],true,true);
}


/*!
  This method evaluates the stabilization term used in immersed boundary
  simulations. According to Mats Larsons suggestion.
*/

bool LinearElasticity::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                                const Vec3& X, const Vec3&) const
{
  if (eKm < 0)
    return true;
  else if (eKm == 0)
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
  int itgType = dualFld.empty() ? STANDARD : ELEMENT_CENTER;
  if (m_mode == SIM::RHS_ONLY) return itgType;

  return itgType | INTERFACE_TERMS | ELEMENT_CORNERS | NORMAL_DERIVS;
}


double LinearElasticity::getThermalStrain (const Vector&, const Vector&,
                                           const Vec3& X) const
{
  if (!myTemp) return 0.0; // No temperature field

  double T0 = myTemp0 ? (*myTemp0)(X) : 0.0;
  double T = (*myTemp)(X);
  return material->getThermalExpansion(T)*(T-T0);
}


bool LinearElasticity::formInitStrainForces (ElmMats& elMat, const Vector& N,
                                             const Matrix& B, const Matrix& C,
                                             const Vec3& X, double detJW) const
{
  if (eS <= 0 || !myTemp)
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


bool LinearElasticity::finalizeElement (LocalIntegral& elmInt,
                                        const FiniteElement& fe,
                                        const TimeDomain& time, size_t iGP)
{
  if (fe.iel > 0)
  {
    size_t iel = fe.iel - 1;
    if (iel < myKmats.size() && eKm > 0)
      myKmats[iel] = static_cast<ElmMats&>(elmInt).A[eKm-1];
    if (iel < myMmats.size() && eM > 0)
      myMmats[iel] = static_cast<ElmMats&>(elmInt).A[eM-1];
  }

  return this->Elasticity::finalizeElement(elmInt,fe,time,iGP);
}
