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
  //! \brief The constructor forwards to the parent class constructor.
  explicit SIMLinElKL(const char* heading, bool shell = false, bool m = false)
    : SIMKLShell(heading,shell), modal(m) {}
  //! \brief Empty destructor.
  virtual ~SIMLinElKL() {}

  //! \brief Returns whether an analytical solution is available or not.
  virtual bool haveAnaSol() const;

protected:
  //! \brief Returns the actual integrand.
  virtual KirchhoffLove* getProblem(int version = 1);

  //! \brief Parses the analytical solution from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parseAnaSol(char* keyWord, std::istream& is);
  //! \brief Parses the analytical solution from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parseAnaSol(const tinyxml2::XMLElement* elem);

  //! \brief Performs some pre-processing tasks on the FE model.
  virtual void preprocessA();

  //! \brief Returns norm index of the integrated volume.
  virtual size_t getVolumeIndex() const;

private:
  bool modal; //!< Modal dynamics simulation flag
};

#endif
