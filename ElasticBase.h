// $Id$
//==============================================================================
//!
//! \file ElasticBase.h
//!
//! \date Jul 04 2014
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class representing FEM integrands for elasticity problems.
//!
//==============================================================================

#ifndef _ELASTIC_BASE_H
#define _ELASTIC_BASE_H

#include "IntegrandBase.h"
#include "Vec3.h"

namespace TimeIntegration { class BDFD2; }


/*!
  \brief Base class representing the FEM integrand of elasticity problems.
*/

class ElasticBase : public IntegrandBase
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  ElasticBase();

public:
  //! \brief The destructor deletes the BDF object, if any.
  virtual ~ElasticBase();

  //! \brief Defines the gravitation vector.
  void setGravity(const Vec3& g) { gravity = g; }
  //! \brief Defines the gravitation vector.
  void setGravity(double gx, double gy = 0.0, double gz = 0.0)
  { gravity.x = gx; gravity.y = gy; gravity.z = gz; }
  //! \brief Returns the gravitation vector.
  const Vec3& getGravity() const { return gravity; }

  //! \brief Defines the number solution vectors.
  void setNoSolutions(size_t n) { nSV = n; }

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Initializes an integration parameter for the integrand.
  //! \param[in] i Index of the integration parameter to define
  //! \param[in] prm The parameter value to assign
  virtual void setIntegrationPrm(unsigned short int i, double prm);
  //! \brief Returns an integration parameter for the integrand.
  //! \param[in] i Index of the integration parameter to return
  virtual double getIntegrationPrm(unsigned short int i) const;

  //! \brief Advances the BDF time step scheme one step forward.
  void advanceStep(double dt, double dtn);

  //! \brief Returns whether this integrand has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return false; }

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld) const { return fld == 1 ? npv : 0; }
  //! \brief Returns the name of a primary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField1Name(size_t i, const char* prefix) const;

  //! \brief Evaluates the point load integrand at a specified point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] pval Load value at the specified point
  virtual bool evalPoint(LocalIntegral& elmInt, const FiniteElement& fe,
                         const Vec3& pval);

  using IntegrandBase::finalizeElement;
  //! \brief Finalizes the element matrices after the numerical integration.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //!
  //! \details This method is used to pass time step size parameters to the
  //! integrand in case of a dynamics simulation, where it is needed to compute
  //! the effective stiffness/mass matrix used in the Newton iterations.
  virtual bool finalizeElement(LocalIntegral& elmInt,
                               const TimeDomain& time, size_t);

protected:
  Vec3 gravity; //!< Gravitation vector

  // Finite element quantities, i.e., indices into element matrices and vectors.
  // These indices will be identical for all elements in a model and can thus
  // be stored here, even when doing multi-threading. Note that the indices are
  // 1-based, since the value zero is used to signal non-existing matrix/vector.
  unsigned short int eKm; //!< Index to element material stiffness matrix
  unsigned short int eKg; //!< Index to element geometric stiffness matrix
  unsigned short int eM;  //!< Index to element mass matrix
  unsigned short int eS;  //!< Index to element load vector
  unsigned short int iS;  //!< Index to element internal force vector
  unsigned short int nSV; //!< Number of consequtive solution vectors in core

  //! \brief Newmark time integration parameters.
  //! \details The interpretation of each parameter
  //! depends on the actual simulator drivers, as follows: <UL>
  //! <LI> 0: Mass-proportional damping coefficient (Rayleigh damping).
  //! <LI> 1: Stiffness-proportional damping coefficient (Rayleigh damping).
  //! <LI> 2: \f$\alpha_H\f$ for nonlinear Newmark drivers.
  //! \f$\beta\f$ or \f$\alpha_m\f$ for linear Newmark drivers. For linear
  //! drivers, a negative value signals that displacement increments
  //! are used as primary unknowns, otherwise accelerations are used.
  //! <LI> 3: A positive value indicates that the solution driver is linear,
  //! and the actual value is then the \f$\gamma\f$ or \f$\alpha_f\f$ parameter.
  //! A zero value indicates a nonlinear driver, with stiffness-proportional
  //! damping (if any) depending on both material and geometric stiffness.
  //! A negative value means a nonlinear driver, with stiffness-proportional
  //! damping (if any), depending on material stiffness only.
  //! <LI> 4: 1.0 if HHTSIM is used, 2.0 if GenAlphaSIM is used, otherwise 0.0.
  double intPrm[5]; //!< </UL>

  TimeIntegration::BDFD2* bdf; //!< BDF time discretization parameters
};

#endif
