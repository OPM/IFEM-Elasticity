// $Id$
//==============================================================================
//!
//! \file ElasticityArgs.C
//!
//! \date May 31 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Preparsing of input files for Elasticity applications.
//!
//==============================================================================

#include "ElasticityArgs.h"
#include "Utilities.h"
#include "tinyxml2.h"


bool ElasticityArgs::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),"eigensolver") &&
      !utl::getAttribute(elem,"mode",eig))
  {
    const char* value = nullptr;
    const tinyxml2::XMLElement* child = elem->FirstChildElement();
    for (; !value && child; child = child->NextSiblingElement())
      value = utl::getValue(child,"mode");
    if (value) eig = atoi(value);
  }

  return this->SIMargsBase::parse(elem);
}
