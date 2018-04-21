// $Id$
//==============================================================================
//!
//! \file SIMShell.C
//!
//! \date Apr 20 2018
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for nonlinear isogeometric FE analysis of shells.
//!
//==============================================================================

#include "SIMShell.h"
#include "NLKirchhoffLoveShell.h"


KirchhoffLove* SIMShell::getProblem (int)
{
  KirchhoffLove* klp = dynamic_cast<KirchhoffLove*>(myProblem);
  if (!myProblem)
    myProblem = klp = new NLKirchhoffLoveShell();

  return klp;
}
