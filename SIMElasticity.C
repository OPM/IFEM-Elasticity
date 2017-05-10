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
#include "SIM3D.h"

//! Plane strain/stress option for 2D problems.
template<> bool SIMElasticity<SIM2D>::planeStrain = false;
//! Axisymmtry option for 2D problems.
template<> bool SIMElasticity<SIM2D>::axiSymmetry = false;
//! Option for output of Gauss points to VTF for 2D problems.
template<> bool SIMElasticity<SIM2D>::GIpointsVTF = false;

//! Dummy plane strain/stress option for 3D problems.
template<> bool SIMElasticity<SIM3D>::planeStrain = false;
//! Dummy axisymmtry option for 3D problems.
template<> bool SIMElasticity<SIM3D>::axiSymmetry = false;
//! Dummy option for output of Gauss points to VTF for 3D problems.
template<> bool SIMElasticity<SIM3D>::GIpointsVTF = false;
