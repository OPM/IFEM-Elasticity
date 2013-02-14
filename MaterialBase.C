// $Id$
//==============================================================================
//!
//! \file MaterialBase.C
//!
//! \date Mar 01 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for material models.
//!
//==============================================================================

#include "MaterialBase.h"
#include "FiniteElement.h"


bool Material::evaluate (Matrix& C, SymmTensor& sigma, double& U,
                         size_t ip, const Vec3& X, const Tensor& F,
                         const SymmTensor& eps, char iop,
                         const TimeDomain* prm, const Tensor* Fpf) const
{
  FiniteElement fe(0,ip);
  return this->evaluate(C,sigma,U,fe,X,F,eps,iop,prm,Fpf);
}
