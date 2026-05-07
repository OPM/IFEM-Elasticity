// $Id$
//==============================================================================
//!
//! \file DruckerPrager.C
//!
//! \date Apr 08 2026
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Linelastic material model with Drucker-Prager yield criterion.
//!
//==============================================================================

#include "DruckerPrager.h"
#include "FiniteElement.h"
#include "Functions.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml2.h"
#include <functional>


DruckerPrager::DruckerPrager(unsigned short int nd, bool ps, bool ax)
  : LinIsotropic(ps,ax), mySigma(nd, ax || !ps)
{
  version = 'C';
  yieldLimit = nullptr;
  alpha = sigy = 0.0;
  Kappa = 1.0;
}


DruckerPrager::~DruckerPrager ()
{
  delete yieldLimit;
}


void DruckerPrager::parse (const tinyxml2::XMLElement* elem)
{
  this->LinIsotropic::parse(elem);

  const tinyxml2::XMLElement* child = elem->FirstChildElement("yieldlimit");
  if (child)
  {
    utl::getAttribute(child,"version",version,false);
    IFEM::cout <<"\n\t  Drucker-Prager version "<< version << std::endl;

    if (const tinyxml2::XMLNode* aval = child->FirstChild(); aval)
    {
      std::string type;
      utl::getAttribute(child,"type",type,true);
      IFEM::cout <<"\t  Yield limit function ";
      yieldLimit = utl::parseTimeFunc(aval->Value(),type);
    }

    if (utl::getAttribute(child,"K",Kappa) && Kappa < 0.778)
    {
      IFEM::cout <<"  ** K = "<< Kappa <<" is outside valid range [0.778,1],";
      Kappa = 0.778;
      IFEM::cout <<" reset to "<< Kappa << std::endl;
    }

    double beta = 35.0;
    utl::getAttribute(child,"beta",beta);
    beta *= M_PI/180.0;
    if (version == 'C')
      alpha = 2.0*sin(beta)/(sqrt(3.0)*(3.0-sin(beta)));
    else
      alpha = tan(beta);
    IFEM::cout <<"\t  alpha = "<< alpha <<"  beta = "<< beta;
  }
}


bool DruckerPrager::evaluate (Matrix& C, SymmTensor& sigma, double& U,
                              const FiniteElement& fe, const Vec3& X,
                              const Tensor& F, const SymmTensor& eps,
                              char iop, const TimeDomain* prm,
                              const Tensor* Fpf) const
{
  bool ok = this->LinIsotropic::evaluate(C,sigma,U,fe,X,F,eps,iop,prm,Fpf);
  if (iop == 1 && ok)
  {
    mySigma = sigma;
    if (Eaging)
      const_cast<DruckerPrager*>(this)->Emod = (*Eaging)(fe.age);
    if (yieldLimit)
      const_cast<DruckerPrager*>(this)->sigy = (*yieldLimit)(fe.age);
  }

  return ok;
}


int DruckerPrager::getNoIntVariables () const
{
  int nvar = 1;
  if (Eaging) nvar += 1;
  if (yieldLimit) nvar += 2;
  return nvar;
}


double DruckerPrager::getInternalVar (int idx, char* label, size_t) const
{
  if (idx > 1 && idx < 3 && !Eaging)
    ++idx;
  if (idx < 1 || idx > this->getNoIntVariables())
  {
    if (label)
      strcpy(label,"zero");
    return 0.0;
  }

  // Lambda function evaluating the Drucker-Prager stress measure, version A.
  std::function<double(const SymmTensor&)> DPstressA =
    [K=Kappa, a=alpha](const SymmTensor& sigma) -> double
  {
    const double sdim = sigma.size() == 4 ? 3.0 : (double)sigma.dim();
    // Calculate the stress invariants (p,q,r)
    // These expressions are taken from the Abaqus Theory Guide
    const double p = -sigma.trace()/sdim; // Hydrostatic pressure
    const SymmTensor S(sigma + p); // Deviatoric stress
    const double q = sqrt(1.5*S.innerProd(S)); // Mises equivalent stress
    if (K >= 1.0) return q - a*p;
    const double r3 = 4.5*S.innerProd(S*S); // 3rd stress invariant
    return q*(0.5+0.5/K-(0.5-0.5/K)*r3/(q*q*q)) - a*p;
  };

  // Lambda function evaluating the Drucker-Prager stress measure, version C.
  std::function<double(const SymmTensor&)> DPstressC =
    [a=alpha](const SymmTensor& sigma) -> double
  {
    const double sdim = sigma.size() == 4 ? 3.0 : (double)sigma.dim();
    // These expressions were given by ChatGPT for uniaxial compression yield.
    const double I1 = sigma.trace(); // First stress invariant
    const SymmTensor S(sigma - I1/sdim); // Deviatoric stress
    const double J2 = 0.5*S.innerProd(S); // Second deviatoric stress invariant
    return a*I1 + sqrt(J2);
  };

  switch (idx) {
  case 1:
    if (label)
      strcpy(label,"sigma(DP)"); // Drucker-Prager stress messure
    else if (version == 'A')
      return DPstressA(mySigma);
    else
      return DPstressC(mySigma);
    break;
  case 2:
    if (label)
      strcpy(label,"Youngs modulus");
    else
      return Emod;
    break;
  case 3:
    if (label)
      strcpy(label,"Yield limit");
    else
      return sigy;
    break;
  case 4:
    if (label)
      strcpy(label,"Yield utilization");
    else if (sigy > 0.0)
    {
      const double osqrt3 = 1.0/sqrt(3.0);
      if (version == 'A')
        return 100.0*DPstressA(mySigma)/((1.0-alpha/3.0)*sigy);
      else
        return 100.0*DPstressC(mySigma)/((alpha+osqrt3)*sigy);
    }
    break;
  }

  return 0.0;
}
