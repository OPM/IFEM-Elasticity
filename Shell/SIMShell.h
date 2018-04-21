// $Id$
//==============================================================================
//!
//! \file SIMShell.h
//!
//! \date Apr 20 2018
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for nonlinear isogeometric FE analysis of shells.
//!
//==============================================================================

#ifndef _SIM_SHELL_H
#define _SIM_SHELL_H

#include "SIMKLShell.h"


/*!
  \brief Driver class for nonlinear isogeometric FE analysis of shells.
*/

class SIMShell : public SIMKLShell
{
public:
  //! \brief Default constructor.
  SIMShell() : SIMKLShell(true) {}
  //! \brief Empty destructor.
  virtual ~SIMShell() {}

protected:
  //! \brief Returns the actual integrand.
  virtual KirchhoffLove* getProblem(int);
};

#endif
