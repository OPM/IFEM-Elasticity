// $Id$
//==============================================================================
//!
//! \file NonlinearElasticityULMixed.h
//!
//! \date Dec 22 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for mixed nonlinear elasticity problems.
//!
//==============================================================================

#ifndef _NONLINEAR_ELASTICITY_UL_MIXED_H
#define _NONLINEAR_ELASTICITY_UL_MIXED_H

#include "NonlinearElasticityUL.h"
#include "ElmMats.h"
#include "Tensor.h"


/*!
  \brief Class representing the integrand of the nonlinear elasticity problem.
  \details This class implements a mixed Updated Lagrangian formulation with
  continuous pressure and volumetric change fields.
*/

class NonlinearElasticityULMixed : public NonlinearElasticityUL
{
  //! \brief Class representing the element matrices of the mixed formulation.
  class MixedElmMats : public ElmMats
  {
  public:
    //! \brief Default constructor.
    MixedElmMats();
    //! \brief Empty destructor.
    virtual ~MixedElmMats() {}
    //! \brief Returns the element-level Newton matrix.
    virtual const Matrix& getNewtonMatrix() const;
    //! \brief Returns the element-level right-hand-side vector
    //! associated with the Newton matrix.
    virtual const Vector& getRHSVector() const;
  };

public:
  //! \brief The default constructor invokes the parent class constructor.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] axS \e If \e true, and axisymmetric 3D formulation is assumed
  NonlinearElasticityULMixed(unsigned short int n = 3, bool axS = false);
  //! \brief Empty destructor.
  virtual ~NonlinearElasticityULMixed() {}

  //! \brief Prints out problem definition to the log stream.
  virtual void printLog() const;

  using NonlinearElasticityUL::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element for each basis
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(const std::vector<size_t>& nen, size_t,
					  bool neumann) const;

  using NonlinearElasticityUL::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Nodal point correspondance for the bases
  //! \param[in] elem_sizes Size of each basis on the element
  //! \param[in] basis_sizes Size of each basis on the patch
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC,
                           const std::vector<size_t>& elem_sizes,
                           const std::vector<size_t>& basis_sizes,
			   LocalIntegral& elmInt);

  using NonlinearElasticityUL::evalIntMx;
  //! \brief Evaluates the mixed field problem integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Mixed finite element data of current integration point
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
			 const TimeDomain& prm, const Vec3& X) const;

  using NonlinearElasticityUL::evalBouMx;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Mixed finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBouMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
			 const Vec3& X, const Vec3& normal) const;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  virtual NormBase* getNormIntegrand(AnaSol* = 0) const;

  //! \brief Returns whether a mixed formulation is used.
  virtual bool mixedFormulation() const { return true; }

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld = 2) const;
  //! \brief Returns the name of a primary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField1Name(size_t i, const char* prefix = 0) const;

  friend class ElasticityNormULMixed;
};


/*!
  \brief Class representing the integrand of the mixed elasticity energy norm.
*/

class ElasticityNormULMixed : public ElasticityNormUL
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The linear elasticity problem to evaluate norms for
  explicit ElasticityNormULMixed(NonlinearElasticityULMixed& p) : ElasticityNormUL(p) {}
  //! \brief Empty destructor.
  virtual ~ElasticityNormULMixed() {}

  using ElasticityNormUL::evalIntMx;
  //! \brief Evaluates the integrand at an interior point (mixed).
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Mixed finite element data of current integration point
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
			 const TimeDomain& prm, const Vec3& X) const;

  using ElasticityNormUL::evalBouMx;
  //! \brief Evaluates the integrand at a boundary point (mixed).
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Mixed finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBouMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
			 const Vec3& X, const Vec3& normal) const;
};

#endif
