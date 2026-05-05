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
#include <map>

class SIMbase;
class SIMgeneric;
class Vec3;


namespace Elastic //! Dimension-independent utilities for Elasticity
{
  //! \brief Prints a norm group to the log stream.
  //! \param[in] gNorm The norm values to print
  //! \param[in] rNorm Reference norms for the first norm group
  //! \param[in] name Projection name associated with this norm group
  //! \param[in] model The FE model associated with the norm values to print
  void printNorms(const Vector& gNorm, const Vector& rNorm,
                  const std::string& name, const SIMbase* model);

  //! \brief Prints out interface force resultants for a set of boundaries.
  //! \param[in] sf Internal nodal forces
  //! \param weights Nodal weights (in case some nodes are present in more sets)
  //! \param[in] bCode Boundary property codes to calculate forces for
  //! \param[in] model The FE model associated with the forces to print
  //! \param[in] indented If \e true, indent the print with two spaces
  void printBoundaryForces(const Vector& sf, RealArray& weights,
                           const std::map<int,Vec3>& bCode,
                           const SIMgeneric* model, bool indented = false);

  extern bool planeStrain; //!< Plane strain/stress option - 2D only
  extern bool axiSymmetry; //!< Axisymmtry option - 2D only
  extern bool GIpointsVTF; //!< Option for Gauss point output to VTF
  extern double time;      //!< Time for function evaluation (linear problems)
}

#endif
