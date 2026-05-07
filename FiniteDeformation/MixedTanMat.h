// $Id$
//==============================================================================
//!
//! \file MixedTanMat.h
//!
//! \date Apr 29 2026
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Utilities for mixed model constitutive tangent matrix.
//!
//==============================================================================

#ifndef _MIXED_TAN_MAT_H
#define _MIXED_TAN_MAT_H

#include "matrix.h"


namespace MixedMat //! Utilities for mixed model constitutive tangent matrix.
{
  //! \brief Extends the (6x6) constitutive matrix by volumetric terms.
  //! \param D Tangent constitutive matrix
  //! \param[in] sigma 3D Cauchy stress tensor
  //! \param[in] p_bar Hydrostatic pressure
  //! \param[in] p_mix Mixed pressure
  void addVolumetricTerms(utl::matrix<double>& D,
                          const std::vector<double>& sigma,
                          double p_bar, double p_mix);
}

#endif
