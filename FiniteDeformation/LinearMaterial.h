// $Id$
//==============================================================================
//!
//! \file LinearMaterial.h
//!
//! \date Mar 08 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief General linear elastic material with push-forward transformations.
//!
//==============================================================================

#ifndef _LINEAR_MATERIAL_H
#define _LINEAR_MATERIAL_H

#include "MaterialBase.h"


/*!
  \brief Class representing a general linear elastic material model.
  \details This class is a wrapper for any linear-elastic material model,
  when it is to be used in a nonlinear finite element solver based on
  the updated Lagrangian formulation. The \a evaluate method of this class
  just invokes the corresponding method of the object pointed to by the
  \a material member, and then transforms the resulting constitutive matrix
  and stress tensor to the current updated reference frame, based on the
  supplied deformation gradient tensor, \a *Fpf (or \a F if \a Fbf is nullptr).
*/

class LinearMaterial : public Material
{
public:
  //! \brief The constructor initializes the material properties pointer.
  //! \param[in] mat Pointer to material model for linear elastic analysis
  explicit LinearMaterial(const Material* mat) : material(mat) {}
  //! \brief The destructor deletes the material properties object.
  virtual ~LinearMaterial() { delete material; }

  //! \brief Returns \e false if plane stress in 2D.
  virtual bool isPlaneStrain() const { return material->isPlaneStrain(); }

  //! \brief Prints out material parameters to the log stream.
  virtual void printLog() const { material->printLog(); }

  //! \brief Evaluates the mass density at an integration point.
  virtual double getMassDensity(const Vec3& X) const;

  using Material::evaluate;
  //! \brief Evaluates the constitutive relation at an integration point.
  //! \param[out] C Constitutive matrix at current point
  //! \param[out] sigma Stress tensor at current point
  //! \param[out] U Strain energy density at current point
  //! \param[in] fe Finite element quantities at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] F Deformation gradient at current point
  //! \param[in] eps Strain tensor at current point
  //! \param[in] iop Calculation option;
  //!  -1 : Calculate the inverse constitutive matrix only,
  //!   0 : Calculate the constitutive matrix only,
  //!   1 : Calculate Cauchy stresses and the constitutive matrix,
  //!   2 : 2nd Piola-Kirchhoff stresses and tangent constitutive matrix,
  //!   3 : Calculate the strain energy density only.
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] Fpf Deformation gradient for push-forward transformation
  virtual bool evaluate(Matrix& C, SymmTensor& sigma, double& U,
                        const FiniteElement& fe, const Vec3& X,
                        const Tensor& F, const SymmTensor& eps,
                        char iop = 1, const TimeDomain* prm = nullptr,
                        const Tensor* Fpf = nullptr) const;

  //! \brief Returns number of internal result variables of the material model.
  virtual int getNoIntVariables() const;

  //! \brief Returns an internal variable associated with the material model.
  virtual double getInternalVar(int idx, char* label, size_t iGP) const;

private:
  const Material* material; //!< Linear-elastic material properties object
};

#endif
