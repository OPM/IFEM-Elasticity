// $Id$
//==============================================================================
//!
//! \file ElasticityArgs.h
//!
//! \date May 31 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Preparsing of input files for Elasticity applications.
//!
//==============================================================================

#ifndef _ELASTICITY_ARGS_H
#define _ELASTICITY_ARGS_H

#include "SIMargsBase.h"


/*!
  \brief Class holding application parameters.
*/

class ElasticityArgs : public SIMargsBase
{
public:
  //! \brief Default constructor.
  ElasticityArgs() : SIMargsBase("elasticity"), eig(0) {}

protected:
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const tinyxml2::XMLElement* elem);

public:
  int eig; //!< Eigenvalue solver option
};

#endif
