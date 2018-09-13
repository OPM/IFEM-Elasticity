// $Id$
//==============================================================================
//!
//! \file IsotropicTextureMat.C
//!
//! \date Sep 13 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Isotropic linear elastic material model defined through a texture.
//!
//==============================================================================

#include "IsotropicTextureMat.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"


void IsotropicTextureMat::parse (const TiXmlElement* elem)
{
  utl::getAttribute(elem, "file", textureFile);
  const TiXmlElement* child = elem->FirstChildElement("range");
  for (; child; child = child->NextSiblingElement("range"))
  {
    std::pair<double,double> minmax;
    utl::getAttribute(child,"min",minmax.first);
    utl::getAttribute(child,"max",minmax.second);
    materials.insert(std::make_pair(minmax,LinIsotropic(planeStress,axiSymmetry)));
    materials[minmax].parse(child);
  }
}


void IsotropicTextureMat::printLog () const
{
  for (const auto it : materials) {
    IFEM::cout << "Material with range [" << it.first.first << "," << it.first.second << "]:\n";
    it.second.printLog();
  }
}


const LinIsotropic* IsotropicTextureMat::findMaterial (const FiniteElement& fe) const
{
  double I = this->findIntensity(fe);
  auto mat = std::find_if(materials.begin(), materials.end(),
            [I](const std::pair<std::pair<double,double>,LinIsotropic>& a)
            {
              return a.first.first >= I && a.first.second <= I;
            });

  if (mat == materials.end())
    return nullptr;

  return &mat->second;
}


double IsotropicTextureMat::findIntensity (const FiniteElement& fe) const
{
  return 0.0;
}


bool IsotropicTextureMat::evaluate (double& lambda, double& mu,
                                    const FiniteElement& fe,
                                    const Vec3& X) const
{
  const LinIsotropic* mat = this->findMaterial(fe);
  if (!mat)
    return false;

  return mat->evaluate(lambda, mu, fe, X);
}


bool IsotropicTextureMat::evaluate (Matrix& C, SymmTensor& sigma, double& U,
                                    const FiniteElement& fe, const Vec3& X,
                                    const Tensor& T, const SymmTensor& eps,
                                    char iop, const TimeDomain* time, const Tensor* T2) const
{
  const LinIsotropic* mat = this->findMaterial(fe);
  if (!mat)
    return false;

  return mat->evaluate(C, sigma, U, fe, X, T, eps, iop, time, T2);
}
