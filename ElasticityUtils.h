// $Id$
//==============================================================================
//!
//! \file ElasticityUtils.h
//!
//! \date Jun 17 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Utilities for Elasticity simulators.
//!
//==============================================================================

#ifndef _ELASTICITY_UTILS_H
#define _ELASTICITY_UTILS_H

#include "MatVec.h"
#include <string>

class SIMbase;


namespace ElasticityUtils //! Dimension-independent utilities for Elasticity
{
  //! \brief Prints a norm group to the log stream.
  //! \param[in] gNorm The norm values to print
  //! \param[in] rNorm Reference norms for the first norm group
  //! \param[in] name Projection name associated with this norm group
  void printNorms(const Vector& gNorm, const Vector& rNorm,
                  const std::string& name, const SIMbase* model);
}

#endif
