// $Id$
//==============================================================================
//!
//! \file LinearElasticity.h
//!
//! \date Nov 12 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for linear elasticity problems.
//!
//==============================================================================

#ifndef _LINEAR_ELASTICITY_H
#define _LINEAR_ELASTICITY_H

#include "Elasticity.h"
#include "GlobalIntegral.h"

class RealFunc;


/*!
  \brief Class representing the integrand of the linear elasticity problem.
  \details Most methods of this class are inherited form the base class.
  Only the \a evalInt method, which is specific for linear elasticity problems
  (and not used in nonlinear problems) is implemented here.
*/

class LinearElasticity : public Elasticity
{
public:
  //! \brief Default constructor.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] axS \e If \e true, an axisymmetric 3D formulation is assumed
  //! \param[in] GPout \e If -e true, write Gauss point coordinates to VTF
  LinearElasticity(unsigned short int n, bool axS = false, bool GPout = false);
  //! \brief Empty destructor.
  virtual ~LinearElasticity() {}

  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  using Elasticity::initIntegration;
  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  //! \param[in] nBp Total number of boundary integration points
  virtual void initIntegration(size_t nGp, size_t nBp);

  using Elasticity::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] X0 Cartesian coordinates of the element center
  //! \param elmInt Local integral for element
  virtual bool initElement(const std::vector<int>& MNPC, const FiniteElement&,
                           const Vec3& X0, size_t, LocalIntegral& elmInt);

  //! \brief Defines the global integral for calculating reaction forces only.
  void setReactionIntegral(GlobalIntegral* q) { delete myReacIt; myReacIt = q; }
  //! \brief Returns the system quantity to be integrated by \a *this.
  virtual GlobalIntegral& getGlobalInt(GlobalIntegral* gq) const;

  //! \brief Returns whether there are any traction values to write to VTF.
  virtual bool hasTractionValues() const;
  //! \brief Writes the surface tractions for a given time step to VTF-file.
  //! \param vtf The VTF-file object to receive the tractions
  //! \param[in] iStep Load/time step identifier
  //! \param geoBlk Running geometry block counter
  //! \param nBlock Running result block counter
  virtual bool writeGlvT(VTF* vtf, int iStep, int& geoBlk, int& nBlock) const;

  using Elasticity::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  //! \brief Evaluates the integrand at an element interface point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Interface normal vector at current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  //! \brief Returns which integrand to be used.
  virtual int getIntegrandType() const;

  //! \brief Returns the initial temperature field.
  const RealFunc* getInitialTemperature() const { return myTemp0; }
  //! \brief Returns the stationary temperature field.
  const RealFunc* getTemperature() const { return myTemp; }

protected:
  //! \brief Evaluates the thermal strain at current integration point.
  //! \param[in] X Cartesian coordinates of current integration point
  virtual double getThermalStrain(const Vector&, const Vector&,
                                  const Vec3& X) const;

  //! \brief Calculates integration point initial strain force contributions.
  //! \param elMat Element matrices for current element
  //! \param[in] N Basis function values at current point
  //! \param[in] B Strain-displacement matrix
  //! \param[in] C Constitutive matrix
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] detJW Jacobian determinant times integration point weight
  virtual bool formInitStrainForces(ElmMats& elMat, const Vector& N,
                                    const Matrix& B, const Matrix& C,
                                    const Vec3& X, double detJW) const;

  RealFunc* myTemp0; //!< Initial temperature field
  RealFunc* myTemp;  //!< Explicit stationary temperature field

private:
  GlobalIntegral*  myReacIt; //!< Reaction-forces-only integral
  mutable Vec3Vec* myItgPts; //!< Global Gauss point coordinates
};

#endif
