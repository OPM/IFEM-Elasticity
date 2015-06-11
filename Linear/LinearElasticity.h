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


/*!
  \brief Class representing the integrand of the linear elasticity problem.
  \details Most methods of this class are inherited form the base class.
  Only the \a evalInt and \a evalBou methods, which are specific for linear
  elasticity problems (and not used in nonlinear problems) are implemented here.
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

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  //! \param[in] nBp Total number of boundary integration points
  virtual void initIntegration(size_t nGp, size_t nBp);

  //! \brief Returns whether there are any traction values to write to VTF.
  virtual bool hasTractionValues() const;
  //! \brief Writes the surface tractions for a given time step to VTF-file.
  //! \param vtf The VTF-file object to receive the tractions
  //! \param[in] iStep Load/time step identifier
  //! \param geoBlk Running geometry block counter
  //! \param nBlock Running result block counter
  virtual bool writeGlvT(VTF* vtf, int iStep, int& geoBlk, int& nBlock) const;

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

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  //! \brief Returns which integrand to be used.
  virtual int getIntegrandType() const;

protected:
  //! \brief Calculates integration point initial strain force contributions.
  virtual bool formInitStrainForces(ElmMats&, const Vector&,
                                    const Matrix&, const Matrix&,
                                    const Vec3&, double) const { return true; }

private:
  mutable Vec3Vec* myItgPts; //!< Global Gauss point coordinates
};

#endif
