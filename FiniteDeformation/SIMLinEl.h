// $Id$
//==============================================================================
//!
//! \file SIMLinEl.h
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#ifndef _SIM_LIN_EL_H
#define _SIM_LIN_EL_H

#include "SIMElasticity.h"
#include "LinearElasticity.h"
#include "ElasticityUtils.h"


/*!
  \brief Driver class for isogeometric FEM analysis of linear elastic problems.
*/

template<class Dim> class SIMLinEl : public SIMElasticity<Dim>
{
public:
  //! \brief Default constructor.
  explicit SIMLinEl(bool checkRHS = false) : SIMElasticity<Dim>(checkRHS) {}
  //! \brief Empty destructor.
  virtual ~SIMLinEl() {}

protected:
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    if (!Dim::myProblem)
      Dim::myProblem = new LinearElasticity(Dim::dimension,
                                            Elastic::axiSymmetry);

    return dynamic_cast<Elasticity*>(Dim::myProblem);
  }
};

#endif
