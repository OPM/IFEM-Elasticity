// $Id$
//==============================================================================
//!
//! \file DruckerPrager.h
//!
//! \date Apr 08 2026
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Linelastic material model with Drucker-Prager yield criterion.
//!
//==============================================================================

#ifndef _DRUCKER_PRAGER_H
#define _DRUCKER_PRAGER_H

#include "LinIsotropic.h"
#include "Tensor.h"


/*!
  \brief Class representing a Drucker-Prager linear elastic material model.
  \details This class extends the LinIsotropic material model with an
  age-dependent yield limit. It also adds the Drucker-Prager stress measure
  as additional result variable together with the age dependent Youngs modulus
  and yield stress limit, as well as the current yield utilization.
*/

class DruckerPrager : public LinIsotropic
{
public:
  //! \brief The default constructor forwards to the parent class constructor.
  //! \param[in] nd Number of spatial dimensions (2 or 3)
  //! \param[in] ps If \e true, assume plane stress in 2D
  //! \param[in] ax If \e true, assume 3D axi-symmetric material
  DruckerPrager(unsigned short int nd = 3, bool ps = false, bool ax = false);
  //! \brief The destructor deletes the yield limit function, if defined.
  virtual ~DruckerPrager();

  //! \brief Parses material parameters from an XML element.
  virtual void parse(const tinyxml2::XMLElement* elem);

  using LinIsotropic::evaluate;
  //! \brief Evaluates the constitutive relation at an integration point.
  //! \param[out] C Constitutive matrix at current point
  //! \param[out] sigma Stress tensor at current point
  //! \param[out] U Strain energy density at current point
  //! \param[in] fe Finite element quantities at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] F Deformation gradient at current point
  //! \param[in] eps Strain tensor at current point
  //! \param[in] iop Calculation option:
  //! - 0 : Calculate the consitutive matrix only
  //! - 1 : Cauchy stresses and the tangent constitutive matrix
  //! - 3 : Calculate strain energy density only
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] Fpf Deformation gradient for push-forward transformation
  virtual bool evaluate(Matrix& C, SymmTensor& sigma, double& U,
                        const FiniteElement& fe, const Vec3& X,
                        const Tensor& F, const SymmTensor& eps, char iop,
                        const TimeDomain* prm, const Tensor* Fpf) const;

  //! \brief Returns number of internal result variables of the material model.
  virtual int getNoIntVariables() const;

  //! \brief Returns an internal variable associated with the material model.
  //! \param[in] idx 1-based index of the internal variable
  //! \param[out] label Name of the internal variable (for result presentation)
  virtual double getInternalVar(int idx, char* label, size_t) const;

private:
  ScalarFunc* yieldLimit; //!< Age-dependent yield limit

  mutable SymmTensor mySigma; //!< Cauchy stress tensor at last evaluation point

  double sigy;  //!< Current yield limit
  double alpha; //!< Material parameter, equals tan(beta)
  double Kappa; //!< Material parameter, valid range [0.778,1.0]
  char version; //!< Formulation version
};

#endif
