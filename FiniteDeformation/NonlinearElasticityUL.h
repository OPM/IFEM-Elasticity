// $Id$
//==============================================================================
//!
//! \file NonlinearElasticityUL.h
//!
//! \date Sep 21 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity problems.
//!
//==============================================================================

#ifndef _NONLINEAR_ELASTICITY_UL_H
#define _NONLINEAR_ELASTICITY_UL_H

#include "Elasticity.h"

class ElmNorm;


/*!
  \brief Class representing the integrand of the nonlinear elasticity problem.
  \details This class implements an Updated Lagrangian formulation. It inherits
  most of the Elasticity methods, but reimplements the \a kinematics method
  for calculating the deformation gradient and the associated Green-Lagrange
  strain tensor. The \a evalInt and \a evalBou methods are also reimplemented
  to account for the updated geometry.
*/

class NonlinearElasticityUL : public Elasticity
{
public:
  //! \brief The constructor invokes the parent class constructor only.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] axS \e If \e true, and axisymmetric 3D formulation is assumed
  //! \param[in] lop Load option (0=on initial length, 1=on updated length)
  NonlinearElasticityUL(unsigned short int n, bool axS = false, char lop = 0);
  //! \brief Empty destructor.
  virtual ~NonlinearElasticityUL() {}

  //! \brief Prints out problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  //! \param[in] nBp Total number of boundary integration points
  virtual void initIntegration(size_t nGp, size_t nBp);
  //! \brief Initializes the integrand for a new integration loop.
  //! \param[in] prm Nonlinear solution algorithm parameters
  virtual void initIntegration(const TimeDomain& prm, const Vector&, bool);
  //! \brief Initializes the integrand for a new result point loop.
  //! \param[in] lambda Load parameter
  //! \param[in] prinDir If &gt; 0, compute/store principal stress directions.
  //! If &lt; 0, Use the integration point buffers for the material instead of
  //! the result point buffers.
  virtual void initResultPoints(double lambda, char prinDir);

  using Elasticity::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const TimeDomain& prm, const Vec3& X) const;

  using Elasticity::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  //!
  //! \details This method is reimplemented in this class to account for
  //! possibly with-rotated traction fields (non-conservative loads).
  //! For uni-directional (conservative) loads, it is similar to the
  //! \a LinearElasticity::evalBou method.
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  virtual NormBase* getNormIntegrand(AnaSol* = nullptr) const;

  //! \brief Returns norm index of the integrated element and total volume.
  virtual size_t getVolumeIndex(bool) const { return 7; }

  //! \brief Calculates some kinematic quantities at current point.
  //! \param[in] eV Element solution vector
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] r Radial coordinate of current point
  //! \param[out] F Deformation gradient at current point
  //! \param[out] E Green-Lagrange strain tensor at current point
  virtual bool kinematics(const Vector& eV,
                          const Vector& N, const Matrix& dNdX, double r,
                          Matrix&, Tensor& F, SymmTensor& E) const;

  //! \brief Returns \e true if simulation diverged on integration point level.
  //! \param[in] iP Global (1-based) index for the integration point to check,
  //! checking all points if \a iP is zero
  virtual bool diverged(size_t iP) const;

private:
  char loadOp; //!< Load option

  friend class ElasticityNormUL;
};


/*!
  \brief Class representing the integrand of the elasticity energy norm.
  \details This class reimplements the \a evalInt method to use the strain
  energy density value returned by the nonlinear constitutive model.
  It also reimplements the \a evalBou method with a path integral of
  the external energy due to boundary tractions.
*/

class ElasticityNormUL : public ElasticityNorm
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The linear elasticity problem to evaluate norms for
  explicit ElasticityNormUL(NonlinearElasticityUL& p) : ElasticityNorm(p) {}
  //! \brief Empty destructor.
  virtual ~ElasticityNormUL() {}

  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  //! \param[in] nBp Total number of boundary integration points
  virtual void initIntegration(size_t nGp, size_t nBp);

  using ElasticityNorm::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const TimeDomain& prm, const Vec3& X) const;

  using ElasticityNorm::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  //! \brief Returns the number of norm quantities.
  virtual size_t getNoFields(int group = 0) const;
  //! \brief Returns the name of a norm quantity.
  std::string getName(size_t i, size_t j, const char* prefix) const;

protected:
  //! \brief Evaluates and accumulates the point-wise norm quantities.
  //! \param pnorm The element norms
  //! \param[in] S Cauchy stress tensor
  //! \param[in] U Strain energy density
  //! \param[in] detF Determinant of deformation gradient
  //! \param[in] detJxW Jacobian determinant times integration point weight
  //!
  //! \details This method is used by the virtual \a evalInt method and is
  //! separated out such that it also can be reused by sub-classes.
  static bool evalInt(ElmNorm& pnorm, const SymmTensor& S,
                      double U, double detF, double detJxW);

private:
  // Data for path-integral of the external energy due to boundary tractions
  mutable RealArray         Ux; //!< External energy densities
  mutable std::vector<Vec3> up; //!< Previous displacements
  mutable std::vector<Vec3> tp; //!< Previous tractions
};

#endif
