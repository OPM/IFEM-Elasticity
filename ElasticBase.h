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

#include "HasGravityBase.h"

class Material;
namespace TimeIntegration { class BDFD2; }


/*!
  \brief Base class representing the FEM integrand of elasticity problems.
*/

class ElasticBase : public HasGravityBase
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  ElasticBase();

public:
  //! \brief The destructor deletes the BDFD2 object, if any.
  virtual ~ElasticBase();

  //! \brief Parses material properties from a character string.
  virtual Material* parseMatProp(char* = nullptr) { return nullptr; }
  //! \brief Parses material properties from an XML-element.
  virtual Material* parseMatProp(const tinyxml2::XMLElement*)
  { return nullptr; }

  //! \brief Defines the material properties.
  virtual void setMaterial(Material*) {}
  //! \brief Returns the current material object.
  virtual Material* getMaterial() const { return nullptr; }

  //! \brief Initializes time integration parameters.
  virtual bool init(const TimeDomain&) { return true; }

  //! \brief Defines the solution mode before the element assembly is started.
  virtual void setMode(SIM::SolutionMode mode);
  //! \brief Returns current solution mode.
  virtual SIM::SolutionMode getMode(bool simMode) const;
  //! \brief Updates the external load vector index for gradient calculation.
  void setLoadGradientMode() { eS = 2; }

  //! \brief Initializes a time integration parameter for the integrand.
  //! \param[in] i Index of the integration parameter to define
  //! \param[in] prm The parameter value to assign
  virtual void setIntegrationPrm(unsigned short int i, double prm);
  //! \brief Returns an integration parameter for the integrand.
  //! \param[in] i Index of the integration parameter to return
  virtual double getIntegrationPrm(unsigned short int i) const;

  //! \brief Returns the total number of solution vectors.
  virtual size_t getNoSolutions(bool allocated = true) const;

  //! \brief Advances the %BDF time step scheme one step forward.
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

  using HasGravityBase::finalizeElement;
  //! \brief Finalizes the element matrices after the numerical integration.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //!
  //! \details This method is used to pass time step size parameters to the
  //! integrand in case of a dynamics simulation, where it is needed to compute
  //! the effective stiffness/mass matrix used in the Newton iterations.
  virtual bool finalizeElement(LocalIntegral& elmInt,
                               const TimeDomain& time, size_t = 0);

protected:
  // Finite element quantities, i.e., indices into element matrices and vectors.
  // These indices will be identical for all elements in a model and can thus
  // be stored here, even when doing multi-threading. Note that the indices are
  // 1-based, since the value zero is used to signal non-existing matrix/vector.
  short int eKm; //!< Index to element material stiffness matrix
  short int eKg; //!< Index to element geometric stiffness matrix
  short int eM;  //!< Index to element mass matrix
  short int eS;  //!< Index to element load vector
  short int gS;  //!< Index to element load gradient vector
  short int iS;  //!< Index to element internal force vector

  unsigned short int nCS; //!< Number of consecutive solution states in core
  unsigned short int nSV; //!< Total number of solution vectors in core

  std::vector<const char*> matNames; //!< Element matrix names (for debug print)
  std::vector<const char*> vecNames; //!< Element vector names (for debug print)

  //! \brief Newmark time integration parameters.
  //! \details The interpretation of each parameter
  //! depends on the actual simulator drivers, as follows: <UL>
  //! <LI>[0]: Mass-proportional damping coefficient (Rayleigh damping).
  //! <LI>[1]: Stiffness-proportional damping coefficient (Rayleigh damping).
  //! <LI>[2]: &alpha;<SUB>H</SUB> for nonlinear Newmark drivers.
  //! &beta; or &alpha;<SUB>m</SUB> for linear Newmark drivers.
  //! For linear drivers, a negative value signals that displacement increments
  //! are used as primary unknowns, otherwise accelerations are used.
  //! <LI>[3]: A positive value indicates that the solution driver is linear,
  //! and the actual value is then the &gamma; or &alpha;<SUB>f</SUB> parameter.
  //! A zero value indicates a nonlinear driver, with stiffness-proportional
  //! damping (if any) depending on both material and geometric stiffness.
  //! A negative value means a nonlinear driver, with stiffness-proportional
  //! damping (if any), depending on material stiffness only.
  //! <LI>[4]: 1.0 if HHTSIM is used, 2.0 if GenAlphaSIM is used, otherwise 0.0.
  double intPrm[5]; //!< </UL>

  TimeIntegration::BDFD2* bdf; //!< BDF time discretization parameters
};

#endif
