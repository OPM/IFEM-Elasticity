// $Id$
//==============================================================================
//!
//! \file NonlinearElasticityULMX.h
//!
//! \date Dec 14 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for nonlinear elasticity mixed problems.
//!
//==============================================================================

#ifndef _NONLINEAR_ELASTICITY_UL_MX_H
#define _NONLINEAR_ELASTICITY_UL_MX_H

#include "NonlinearElasticityUL.h"


/*!
  \brief Class representing the integrand of the nonlinear elasticity problem.
  \details This class implements a mixed Updated Lagrangian formulation with
  internal pressure and volumetric change modes.
*/

class NonlinearElasticityULMX : public NonlinearElasticityUL
{
public:
  //! \brief The default constructor invokes the parent class constructor.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] axS \e If \e true, and axisymmetric 3D formulation is assumed
  //! \param[in] pp Polynomial order of the pressure/volumetric-change field
  NonlinearElasticityULMX(unsigned short int n = 3,
			  bool axS = false, int pp = 1);
  //! \brief Empty destructor.
  virtual ~NonlinearElasticityULMX() {}

  //! \brief Prints out problem definition to the log stream.
  virtual void printLog() const;

  using NonlinearElasticityUL::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
					  bool neumann) const;

  using NonlinearElasticityUL::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] Xc Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points in this element
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC, const FiniteElement&,
			   const Vec3& Xc, size_t nPt,
			   LocalIntegral& elmInt);

  using NonlinearElasticityUL::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
		       const TimeDomain& prm, const Vec3& X) const;

  using NonlinearElasticityUL::finalizeElement;
  //! \brief Finalizes the element matrices after the numerical integration.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] iG Global index of the first integration point in the element
  virtual bool finalizeElement(LocalIntegral& elmInt,
                               const TimeDomain& prm, size_t iG);

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  virtual NormBase* getNormIntegrand(AnaSol* = 0) const;

  //! \brief Returns which integrand to be used.
  virtual int getIntegrandType() const { return ELEMENT_CENTER; }
  //! \brief Returns the number of reduced-order integration points.
  //! \return -1 which is used to flag that the integrand only needs to know
  //! the number of integration points within each element.
  //! \sa The argument \a nPt of the initElement() method.
  virtual int getReducedIntegration(int) const { return -1; }

private:
  int p; //!< Polynomial order of the internal pressure field

  friend class ElasticityNormULMX;
};


/*!
  \brief Class representing the integrand of the mixed elasticity energy norm.
*/

class ElasticityNormULMX : public ElasticityNormUL
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The linear elasticity problem to evaluate norms for
  explicit ElasticityNormULMX(NonlinearElasticityULMX& p) : ElasticityNormUL(p){}
  //! \brief Empty destructor.
  virtual ~ElasticityNormULMX() {}

  using ElasticityNormUL::getLocalIntegral;
  //! \brief Returns a local integral contribution object for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] iEl The element number
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t iEl,
					  bool neumann) const;

  using ElasticityNormUL::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] Xc Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points in this element
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC, const FiniteElement&,
			   const Vec3& Xc, size_t nPt,
			   LocalIntegral& elmInt);

  using ElasticityNormUL::initElementBou;
  //! \brief Initializes current element for boundary integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt The local integral object for current element
  virtual bool initElementBou(const std::vector<int>& MNPC,
			      LocalIntegral& elmInt);

  using ElasticityNormUL::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
		       const TimeDomain& prm, const Vec3& X) const;

  using ElasticityNormUL::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
		       const Vec3& X, const Vec3& normal) const;

  using ElasticityNormUL::finalizeElement;
  //! \brief Finalizes the element norms after the numerical integration.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] iG Global index of the first integration point in the element
  virtual bool finalizeElement(LocalIntegral& elmInt,
                               const TimeDomain& prm, size_t iG);
};

#endif
