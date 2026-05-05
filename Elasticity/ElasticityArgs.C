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


bool ElasticityArgs::parseArg (const char* argv)
{
  if (this->SIMargsBase::parseArg(argv))
    return true;
  else if (!strcmp(argv,"-sup"))
    dim = 4;
  else if (!strcmp(argv,"-1D3D"))
    dim = 13;
  else if (!strcmp(argv,"-1Dsup"))
    dim = 14;
  else
    return false;

  return true;
}


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
  else if (!strcasecmp(elem->Value(),"superelem") && dim == 3)
    dim = 4;
  else if (!strcasecmp(elem->Value(),"beam") && dim == 4)
    dim = 14;

  return this->SIMargsBase::parse(elem);
}
