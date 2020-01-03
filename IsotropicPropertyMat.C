// $Id$
//==============================================================================
//!
//! \file IsotropicPropertyMat.C
//!
//! \date Jan 3 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Isotropic linear elastic material model defined through a table.
//!
//==============================================================================

#include "IsotropicPropertyMat.h"
#include "FiniteElement.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"
#include "StbImage.h"


void IsotropicPropertyMat::parse (const TiXmlElement* elem)
{
  this->TextureProperties::parse(elem);
  if (this->hasProperty("stiffness"))
    Efunc = new PropertyFunc("stiffness", *this);
  if (this->hasProperty("poisson"))
    nuF = new PropertyFunc("poisson", *this);
}


void IsotropicPropertyMat::printLog () const
{
  IFEM::cout << "Isotropic material based on property tables";
  this->TextureProperties::printLog();
}
