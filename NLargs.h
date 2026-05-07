// $Id$
//==============================================================================
//!
//! \file NLargs.h
//!
//! \date Nov 20 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Pre-parsing of input files for Finite Deformation applications.
//!
//==============================================================================

#ifndef _NL_ARGS_H
#define _NL_ARGS_H

#include "XMLInputBase.h"


//! \brief Enum defining various solution drivers.
enum Algorithm { ARCLEN, STATIC, GENALPHA, OLDHHT, NEWHHT };


/*!
  \brief Class for input file pre-parsing and command-line argument values.
*/

class NLargs : public XMLInputBase
{
public:
  //! \brief Parses a command-line argument.
  bool parseArg(const char* argv);

protected:
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const tinyxml2::XMLElement* elem);

private:
  //! \brief Parses the formulation tag from an XML element.
  void parseFormulation(const tinyxml2::XMLElement* elem);
  //! \brief Updates the problem formulation options.
  void setFormulation(int form, int pOrd);

public:
  std::vector<int> options; //!< Problem formulation options
  Algorithm algor = STATIC; //!< Solution algorithm driver flag
  bool adaptive = false; //!< Run and adaptive simulation?
  bool checkRHS = false; //!< Check for right-hand-side patch orientation?
  bool fixDup   = false; //!< Check and resolve duplicated nodes?
  bool twoD     = false; //!< Use the 2D driver?
  char printMax = false; //!< Print out max field values?
  bool dNodeMap = false; //!< Dump node mapping to HDF5 file?
};

#endif
