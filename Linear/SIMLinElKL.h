// $Id$
//==============================================================================
//!
//! \file SIMLinElKL.h
//!
//! \date Sep 16 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based FE analysis of Kirchhoff-Love shells.
//!
//==============================================================================

#ifndef _SIM_LIN_EL_KL_H
#define _SIM_LIN_EL_KL_H

#include "SIMKLShell.h"


/*!
  \brief Driver class for isogeometric FEM analysis of Kirchhoff-Love shells.
*/

class SIMLinElKL : public SIMKLShell
{
public:
  //! \brief Default constructor.
  SIMLinElKL(bool shell = false) : SIMKLShell(shell) {}
  //! \brief Empty destructor.
  virtual ~SIMLinElKL() {}

  //! \brief Returns whether an analytical solution is available or not.
  virtual bool haveAnaSol() const;

protected:
  //! \brief Returns the actual integrand.
  virtual KirchhoffLove* getProblem(int version = 1);

  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Performs some pre-processing tasks on the FE model.
  virtual void preprocessA();
};

#endif
