// $Id$
//==============================================================================
//!
//! \file NonlinearElasticityTL.h
//!
//! \date May 25 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity problems.
//!
//==============================================================================

#ifndef _NONLINEAR_ELASTICITY_TL_H
#define _NONLINEAR_ELASTICITY_TL_H

#include "Elasticity.h"


/*!
  \brief Class representing the integrand of the nonlinear elasticity problem.
  \details This class implements a Total Lagrangian formulation in matrix form.
  It inherits most of the Elasticity methods, but reimplements the \a kinematics
  method, for calculating the nonlinear variant of the strain-displacement
  matrix, \b B, and the associated Green-Lagrange strain tensor, \b E.
  The \a evalBou method is also reimplemented to account for with-rotated loads.
*/

class NonlinearElasticityTL : public Elasticity
{
public:
  //! \brief The constructor invokes the parent class constructor only.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] axS If \e true, an axisymmetric 3D formulation is assumed
  //! \param[in] evalAtElmCenter If \e true, evaluate material at element center
  NonlinearElasticityTL(unsigned short int n, bool axS = false,
                        bool evalAtElmCenter = false);
  //! \brief Empty destructor.
  virtual ~NonlinearElasticityTL() {}

  //! \brief Prints out problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const { return myIntegrandType; }

  using Elasticity::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  using Elasticity::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  //!
  //! \details This method is reimplemented in this class to account for
  //! possibly with-rotated traction fields in the Total-Lagrangian setting.
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

protected:
  //! \brief Calculates some kinematic quantities at current point.
  //! \param[in] eV Element solution vector
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] r Radial coordinate of current point
  //! \param[out] F Deformation gradient at current point
  //! \param[out] Bmat The strain-displacement matrix
  //! \param[out] E Green-Lagrange strain tensor at current point
  //!
  //! \details The deformation gradient \b F and the nonlinear
  //! strain-displacement matrix \b B are established.
  //! The B-matrix is formed only when the variable \a formB is true.
  virtual bool kinematics(const Vector& eV,
                          const Vector& N, const Matrix& dNdX, double r,
                          Matrix& Bmat, Tensor& F, SymmTensor& E) const;

private:
  Integrand::Traits myIntegrandType; //!< Defines additional terms to be used

protected:
  bool formB; //!< Flag determining whether we need to form the B-matrix
};

#endif
