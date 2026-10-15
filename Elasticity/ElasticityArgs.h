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
  \brief Class holding Elasticity application parameters.
*/

class ElasticityArgs : public SIMargsBase
{
public:
  //! \brief Default constructor.
  ElasticityArgs() : SIMargsBase("elasticity") {}

  //! \brief Parses a command-line argument.
  virtual bool parseArg(const char* argv);

protected:
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const tinyxml2::XMLElement* elem);

public:
  int  eig = 0;     //!< Eigenvalue solver option
  bool KLp = false; //!< If \e true, this is  Kirchhoff-Love plate/shell problem
  bool isC1 = false; //!< If \e true, a C1-continuous descretization is used
  bool shell = false; //!< If \e true, this is a shell problem
};

#endif
