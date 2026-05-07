// $Id$
//==============================================================================
//!
//! \file NLoptions.C
//!
//! \date Mar 22 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class for encapsulation of nonlinear formulation options.
//!
//==============================================================================

#include "NLoptions.h"
#include "NonlinearElasticityFbar.h"
#include "NonlinearElasticityULMixed.h"
#include "NonlinearElasticityULMX.h"
#include "NonlinearElasticity.h"
#include "LinearElasticity.h"
#include "ElasticityUtils.h"


ElasticBase* NLoptions::getIntegrand () const
{
  switch (form)
    {
    case SIM::FBAR:
      // F-bar formulation
      return new NonlinearElasticityFbar(ndim,Elastic::axiSymmetry,pOrd);

    case SIM::MIXED_QnQn1:
      // Continuous volumetric change and pressure fields
      return new NonlinearElasticityULMixed(ndim,Elastic::axiSymmetry);

    case SIM::MIXED_QnPn1:
      // Local discontinuous volumetric change and pressure fields
      return new NonlinearElasticityULMX(ndim,Elastic::axiSymmetry,pOrd);

    case SIM::UPDATED_LAGRANGE:
      return new NonlinearElasticityUL(ndim,Elastic::axiSymmetry);

    case SIM::TOTAL_LAGRANGE:
      return new NonlinearElasticityTL(ndim,Elastic::axiSymmetry);

    case SIM::NONLINEAR:
      // Old tensor-based TL-formulation
      return new NonlinearElasticity(ndim);

    case SIM::LINEAR:
      // Small-strain linear formulation
      return new LinearElasticity(ndim,Elastic::axiSymmetry,false,true);

    default:
      std::cerr <<" *** NLoptions::getIntegrand: Unknown problem formulation ("
                << form <<")."<< std::endl;
    }

  return nullptr;
}
