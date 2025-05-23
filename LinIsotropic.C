// $Id$
//==============================================================================
//!
//! \file LinIsotropic.C
//!
//! \date Mar 01 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Isotropic linear elastic material model.
//!
//==============================================================================

#include "LinIsotropic.h"
#include "FiniteElement.h"
#include "Field.h"
#include "Functions.h"
#include "Utilities.h"
#include "Tensor.h"
#include "Vec3.h"
#include "IFEM.h"
#include "tinyxml2.h"


LinIsotropic::LinIsotropic (bool ps, bool ax) : planeStress(ps), axiSymmetry(ax)
{
  Efunc  = nuFunc = rhoFunc = nullptr;
  Efield = nullptr;
  Cpfunc = Afunc = condFunc = nullptr;

  // Default material properties - typical values for steel (SI units)
  Emod = 2.05e11;
  nu = 0.29;
  rho = 7.85e3;
  alpha = 1.2e-7;
  heatcapacity = conductivity = 1.0;
}


LinIsotropic::LinIsotropic (RealFunc* E, double v, double den, bool ps, bool ax)
  : Efunc(E), Efield(nullptr), nu(v), rho(den), planeStress(ps), axiSymmetry(ax)
{
  nuFunc = rhoFunc = nullptr;
  Cpfunc = Afunc = condFunc = nullptr;

  Emod = -1.0; // Should not be referenced
  alpha = 1.2e-7;
  heatcapacity = conductivity = 1.0;
}


LinIsotropic::LinIsotropic (Field* E, double v, double den, bool ps, bool ax)
  : Efunc(nullptr), Efield(E), nu(v), rho(den), planeStress(ps), axiSymmetry(ax)
{
  nuFunc = rhoFunc = nullptr;
  Cpfunc = Afunc = condFunc = nullptr;

  Emod = -1.0; // Should not be referenced
  alpha = 1.2e-7;
  heatcapacity = conductivity = 1.0;
}


LinIsotropic::~LinIsotropic ()
{
  delete Efield;
  delete Efunc;
  delete nuFunc;
  delete rhoFunc;
  delete Afunc;
  delete Cpfunc;
  delete condFunc;
}


void LinIsotropic::parse (const tinyxml2::XMLElement* elem)
{
  if (Emod >= 0.0 && utl::getAttribute(elem,"E",Emod))
    IFEM::cout <<" "<< Emod;
  if (utl::getAttribute(elem,"nu",nu))
    IFEM::cout <<" "<< nu;
  if (utl::getAttribute(elem,"rho",rho))
    IFEM::cout <<" "<< rho;
  if (utl::getAttribute(elem,"alpha",alpha))
    IFEM::cout <<" "<< alpha;
  if (utl::getAttribute(elem,"cp",heatcapacity))
    IFEM::cout <<" "<< heatcapacity;
  if (utl::getAttribute(elem,"kappa",conductivity))
    IFEM::cout <<" "<< conductivity;

  // Lambda function for parsing a spatial property function.
  auto&& parseSpatialFunc = [](const tinyxml2::XMLElement* child, const char* name)
  {
    std::string type;
    utl::getAttribute(child,"type",type,true);
    IFEM::cout <<"\n\t  "<< name <<" function ("<< type <<") ";
    const tinyxml2::XMLNode* aval = child->FirstChild();
    return aval ? utl::parseRealFunc(aval->Value(),type) : nullptr;
  };

  // Lambda function for parsing a scalar property function.
  auto&& parseScalarFunc = [](const tinyxml2::XMLElement* child)
  {
    std::string type;
    utl::getAttribute(child,"type",type,true);
    IFEM::cout <<" ";
    const tinyxml2::XMLNode* aval = child->FirstChild();
    return aval ? utl::parseTimeFunc(aval->Value(),type) : nullptr;
  };

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (Emod >= 0.0 && !Efunc && !strcasecmp(child->Value(),"stiffness"))
      Efunc = parseSpatialFunc(child,"Stiffness");
    else if (!strcasecmp(child->Value(),"poisson"))
      nuFunc = parseSpatialFunc(child,"Poisson's ratio");
    else if (!strcasecmp(child->Value(),"density"))
      rhoFunc = parseSpatialFunc(child,"Mass density");
    else if (!strcasecmp(child->Value(),"thermalexpansion"))
      Afunc = parseScalarFunc(child);
    else if (!strcasecmp(child->Value(),"heatcapacity"))
      Cpfunc = parseScalarFunc(child);
    else if (!strcasecmp(child->Value(),"conductivity"))
      condFunc = parseScalarFunc(child);

  if (!Efunc && !nuFunc && !rhoFunc && !Afunc && !Cpfunc && !condFunc)
    IFEM::cout << std::endl;
}


void LinIsotropic::printLog () const
{
  IFEM::cout <<"LinIsotropic: ";
  if (axiSymmetry)
    IFEM::cout <<"axial-symmetric, ";
  else if (planeStress)
    IFEM::cout <<"plane stress, ";
  IFEM::cout <<"E = "<< Emod <<", nu = "<< nu <<", rho = "<< rho
             <<", alpha = "<< alpha << std::endl;
}


/*!
  The consitutive matrix for Isotropic linear elastic problems
  is defined as follows:

  For 2D plain stress: \f[
  [C] = \frac{E}{(1-\nu^2)} \left[\begin{array}{ccc}
  1 & \ \ \nu & 0 \\
  \nu & \ \ 1 & 0 \\
  0 & \ \ 0 & \frac{1}{2}(1-\nu)
  \end{array}\right] \f]

  For 2D plain strain: \f[
  [C] = \frac{E}{(1+\nu)(1-2\nu)} \left[\begin{array}{ccc}
  1-\nu & \nu & 0 \\
  \nu & 1-\nu & 0 \\
  0 & 0 & \frac{1}{2}-\nu
  \end{array}\right] \f]

  For 3D axisymmetric solids: \f[
  [C] = \frac{E}{(1+\nu)(1-2\nu)} \left[\begin{array}{cccc}
  1-\nu & \nu & \nu & 0 \\
  \nu & 1-\nu & \nu & 0 \\
  \nu & \nu & 1-\nu & 0 \\
  0 & 0 & 0 & \frac{1}{2}-\nu
  \end{array}\right] \f]

  For 3D: \f[
  [C] = \frac{E}{(1+\nu)(1-2\nu)} \left[\begin{array}{cccccc}
  1-\nu & \nu & \nu & 0 & 0 & 0 \\
  \nu & 1-\nu & \nu & 0 & 0 & 0 \\
  \nu & \nu & 1-\nu & 0 & 0 & 0 \\
  0 & 0 & 0 & \frac{1}{2}-\nu & 0 & 0 \\
  0 & 0 & 0 & 0 & \frac{1}{2}-\nu & 0 \\
  0 & 0 & 0 & 0 & 0 & \frac{1}{2}-\nu
  \end{array}\right] \f]
*/

bool LinIsotropic::evaluate (Matrix& C, SymmTensor& sigma, double& U,
                             const FiniteElement& fe, const Vec3& X,
                             const Tensor&, const SymmTensor& eps, char iop,
                             const TimeDomain*, const Tensor*) const
{
  const size_t nsd = sigma.dim();
  const size_t nst = nsd == 2 && axiSymmetry ? 4 : nsd*(nsd+1)/2;
  C.resize(nst,nst,true);

  // Evaluate the scalar stiffness function or field, if defined
  double E = Emod;
  if (Efield)
    E = Efield->valueFE(fe);
  else if (Efunc)
    E = (*Efunc)(X);

  if (nuFunc)
    const_cast<LinIsotropic*>(this)->nu = (*nuFunc)(X);

  if (nsd == 1)
  {
    // Special for 1D problems
    C(1,1) = iop < 0 ? 1.0/E : E;
    if (iop > 0)
    {
      sigma = eps; sigma *= E;
      if (iop == 3)
        U = 0.5*sigma(1,1)*eps(1,1);
    }
    return true;
  }
  else if (nu < 0.0 || nu >= 0.5)
  {
    std::cerr <<" *** LinIsotropic::evaluate: Poisson's ratio "<< nu
              <<" out of range [0,0.5>."<< std::endl;
    return false;
  }

  if (iop < 0) // The inverse C-matrix is wanted
    if (nsd == 3 || (nsd == 2 && (planeStress || axiSymmetry)))
    {
      C(1,1) = 1.0 / E;
      C(2,1) = -nu / E;
    }
    else // 2D plain strain
    {
      C(1,1) = (1.0 - nu*nu) / E;
      C(2,1) = (-nu - nu*nu) / E;
    }

  else
    if (nsd == 2 && planeStress && !axiSymmetry)
    {
      C(1,1) = E / (1.0 - nu*nu);
      C(2,1) = C(1,1) * nu;
    }
    else // 2D plain strain, axisymmetric or 3D
    {
      double fact = E / ((1.0 + nu) * (1.0 - nu - nu));
      C(1,1) = fact * (1.0 - nu);
      C(2,1) = fact * nu;
    }

  C(1,2) = C(2,1);
  C(2,2) = C(1,1);

  const double G = E / (2.0 + nu + nu);
  C(nsd+1,nsd+1) = iop < 0 ? 1.0 / G : G;

  if (nsd == 2 && axiSymmetry)
  {
    C(4,4) = C(3,3);
    C(3,1) = C(2,1);
    C(3,2) = C(2,1);
    C(1,3) = C(2,1);
    C(2,3) = C(2,1);
    C(3,3) = C(1,1);
  }
  else if (nsd > 2)
  {
    C(3,1) = C(2,1);
    C(3,2) = C(2,1);
    C(1,3) = C(2,1);
    C(2,3) = C(2,1);
    C(3,3) = C(1,1);
    C(5,5) = C(4,4);
    C(6,6) = C(4,4);
  }

  if (iop > 0)
  {
    // Calculate the stress tensor, sigma = C*eps
    RealArray sig; // Use a local variable to avoid redimensioning of sigma
    if (eps.dim() != sigma.dim())
    {
      // Account for non-matching tensor dimensions
      SymmTensor epsil(sigma.dim(), nsd == 2 && axiSymmetry);
      if (!C.multiply(epsil=eps,sig))
        return false;
    }
    else
      if (!C.multiply(eps,sig))
        return false;

    sigma = sig; // Add sigma_zz in case of plane strain
    if (!planeStress && ! axiSymmetry && nsd == 2 && sigma.size() == 4)
      sigma(3,3) = nu * (sigma(1,1)+sigma(2,2));
  }

  if (iop == 3) // Calculate strain energy density, // U = 0.5*sigma:eps
    U = 0.5*sigma.innerProd(eps);

  return true;
}


bool LinIsotropic::evaluate (double& lambda, double& mu,
                             const FiniteElement& fe, const Vec3& X) const
{
  if (nuFunc)
    const_cast<LinIsotropic*>(this)->nu = (*nuFunc)(X);

  if (nu < 0.0 || nu >= 0.5)
  {
    std::cerr <<" *** LinIsotropic::evaluate: Poisson's ratio "<< nu
              <<" out of range [0,0.5>."<< std::endl;
    return false;
  }

  // Evaluate the scalar stiffness function or field, if defined
  double E = Emod;
  if (Efield)
    E = Efield->valueFE(fe);
  else if (Efunc)
    E = (*Efunc)(X);

  // Evaluate the Lame parameters
  mu = 0.5*E/(1.0+nu);
  lambda = mu*nu/(0.5-nu);

  return true;
}


double LinIsotropic::getStiffness (const Vec3& X) const
{
  return Efunc ? (*Efunc)(X) : Emod;
}


double LinIsotropic::getPlateStiffness (const Vec3& X, double t) const
{
  double E = Efunc ? (*Efunc)(X) : Emod;
  double v = nuFunc ? (*nuFunc)(X) : nu;
  return E*t*t*t / (12.0 - 12.0*v*v);
}


double LinIsotropic::getMassDensity (const Vec3& X) const
{
  return rhoFunc ? (*rhoFunc)(X) : rho;
}


double LinIsotropic::getThermalExpansion (double T) const
{
  return Afunc ? (*Afunc)(T) : alpha;
}


double LinIsotropic::getHeatCapacity (double T) const
{
  return Cpfunc ? (*Cpfunc)(T) : heatcapacity;
}


double LinIsotropic::getThermalConductivity (double T) const
{
  return condFunc ? (*condFunc)(T) : conductivity;
}
