// $Id$
//==============================================================================
//!
//! \file SIMElasticity.C
//!
//! \date Dec 04 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#include "SIMElasticity.h"
#include "SIM2D.h"

//! Plane strain/stress option for 2D problems.
template<> bool SIMElasticity<SIM2D>::planeStrain = false;
//! Axisymmtry option for 2D problems.
template<> bool SIMElasticity<SIM2D>::axiSymmetry = false;
//! Option for output of Gauss points to VTF for 2D problems.
template<> bool SIMElasticity<SIM2D>::GIpointsVTF = false;
