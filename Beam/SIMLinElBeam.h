// $Id$
//==============================================================================
//!
//! \file SIMLinElBeam.h
//!
//! \date Mar 17 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution drivers for linear elastic beam solvers.
//!
//==============================================================================

#ifndef _SIM_LIN_EL_BEAM_H
#define _SIM_LIN_EL_BEAM_H

#ifdef HAS_GEOMETRY
#include "SIMBeamGeometry.h"
using SIMLinElBeam = SIMBeamGeometry; //!< Beam model with tesselated geometry
#else
#include "SIMElasticBar.h"
using SIMLinElBeam = SIMElasticBar; //!< Standard beam model
#endif

#endif
