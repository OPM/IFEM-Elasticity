// $Id$
//==============================================================================
//!
//! \file SIMElasticity.h
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#ifndef _SIM_ELASTICITY_INSTANCE_H
#define _SIM_ELASTICITY_INSTANCE_H

#include "SIMElasticity.h"
#include "SIM2D.h"
#include "SIM3D.h"

template<> bool SIMElasticity<SIM2D>::planeStrain;
template<> bool SIMElasticity<SIM3D>::planeStrain;
template<> bool SIMElasticity<SIM2D>::axiSymmetry;
template<> bool SIMElasticity<SIM3D>::axiSymmetry;

#endif
