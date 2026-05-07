// $Id$
//==============================================================================
//!
//! \file NeoHookeMaterial.h
//!
//! \date Mar 08 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Neo-Hookean hyperelastic material model.
//!
//==============================================================================

#ifndef _NEO_HOOKE_MATERIAL_H
#define _NEO_HOOKE_MATERIAL_H

#include "LinIsotropic.h"


/*!
  \brief Class representing a Neo-Hookean hyperelastic material model.

  \details This class implements two versions of the Neo-Hookean material model.
  Which version to use is governed by the \ref mVER parameter.
  Only isotropic material properties are supported by this class.
  In 2D, only plane strain is supported.
*/

class NeoHookeMaterial : public LinIsotropic
{
public:
  //! \brief Default constructor.
  explicit NeoHookeMaterial(int ver = 1);
  //! \brief Constructor initializing the material parameters.
  NeoHookeMaterial(double E, double v, double density, int ver = 1);

  //! \brief Parses material parameters from an XML element.
  virtual void parse(const tinyxml2::XMLElement* elem);

  //! \brief Prints out material parameters to the log stream.
  virtual void printLog() const;

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
  //! - 0 : Calculate the constitutive matrix only
  //! - 1 : Cauchy stresses and the tangent constitutive matrix
  //! - 2 : 2nd Piola-Kirchhoff stresses and the tangent constitutive matrix
  //! - 3 : Calculate strain energy density only
  //! \param[in] Fpf Deformation gradient for push-forward transformation
  virtual bool evaluate(Matrix& C, SymmTensor& sigma, double& U,
                        const FiniteElement& fe, const Vec3& X,
                        const Tensor& F, const SymmTensor& eps,
                        char iop = 1, const TimeDomain* = nullptr,
                        const Tensor* Fpf = nullptr) const;

  //! \brief Returns number of internal result variables of the material model.
  virtual int getNoIntVariables() const { return 1; }

  //! \brief Returns an internal variable associated with the material model.
  //! \param[out] label Name of the internal variable (for result presentation)
  virtual double getInternalVar(int, char* label, size_t) const;

protected:
  //! \brief Calculates Lame parameters from \a E and &nu;, or vice versa.
  void findLameParams(double E = -1.0);

  //! \brief Performs calculations for the standard NeoHookean material model.
  //! \param[in] J Determinant of deformation gradient
  //! \param S Left Cauchy-Green deformation tensor / Cauchy stress tensor
  //! \param[out] C Constitutive tensor
  //! \return Strain energy density
  double stdNeoHooke(double J, SymmTensor& S, Matrix& C) const;

  //! \brief Performs calculations for the modified NeoHookean material model.
  //! \param[in] J Determinant of deformation gradient
  //! \param S Left Cauchy-Green deformation tensor / Cauchy stress tensor
  //! \param[out] C Constitutive tensor
  //! \return Strain energy density
  double modNeoHooke(double J, SymmTensor& S, Matrix& C) const;

  //! \brief Computes the volumetric material moduli.
  //! \param[in] J Determinant of deformation gradient
  //! \param[out] U Strain energy density
  //! \param[out] Up Mean stress (pressure)
  //! \param[out] Upp Volumetric constitutive tensor
  void volumetricModuli(double J, double& U, double& Up, double& Upp) const;

private:
  int    mVER; //!< Material version
  int    iVOL; //!< Volumetric option
  double Bmod; //!< Bulk modulus (Lame parameter &lambda;)
  double Smod; //!< Shear modulus (Lame parameter &mu;)

  mutable double sigma_p; //!< Hydrostatic pressure at last evaluation point
};

#endif
