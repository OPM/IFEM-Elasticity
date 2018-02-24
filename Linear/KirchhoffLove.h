// $Id$
//==============================================================================
//!
//! \file KirchhoffLove.h
//!
//! \date Sep 13 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for linear Kirchhoff-Love thin plate and shell problems.
//!
//==============================================================================

#ifndef _KIRCHHOFF_LOVE_H
#define _KIRCHHOFF_LOVE_H

#include "IntegrandBase.h"
#include "Vec3.h"

class LocalSystem;
class Material;
class RealFunc;


/*!
  \brief Class representing the integrand of thin plate and shell problems.

  \details This class contains some problem parameters and common methods
  for representing a linear thin plate or shell problem, based on Kirchhoff-Love
  theory. The actual integrand implementation is located in sub-classes.
*/

class KirchhoffLove : public IntegrandBase
{
public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions (1=beam, 2=plate, 3=shell)
  explicit KirchhoffLove(unsigned short int n = 2);
  //! \brief The destructor frees the dynamically allocated data objects.
  virtual ~KirchhoffLove();

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Defines the gravitation constant.
  void setGravity(double g) { gravity = g; }
  //! \brief Defines the plate/shell thickness.
  void setThickness(double t) { thickness = t; }
  //! \brief Defines the pressure field.
  void setPressure(RealFunc* pf) { presFld = pf; }
  //! \brief Defines the material properties.
  void setMaterial(Material* mat) { material = mat; }

  //! \brief Defines the local coordinate system for stress resultant output.
  void setLocalSystem(LocalSystem* cs) { locSys = cs; }

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const { return SECOND_DERIVATIVES; }

  using IntegrandBase::initIntegration;
  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  //! \param[in] nBp Total number of boundary integration points
  virtual void initIntegration(size_t nGp, size_t nBp);

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

  //! \brief Evaluates the pressure field (if any) at specified point.
  virtual double getPressure(const Vec3& X) const;
  //! \brief Returns whether an external load is defined.
  virtual bool haveLoads() const;

  //! \brief Returns the derivative order of the differential operator.
  virtual int derivativeOrder() const { return 2; }

  //! \brief Writes the surface pressure for a given time step to VTF-file.
  //! \param vtf The VTF-file object to receive the pressure vectors
  //! \param[in] iStep Load/time step identifier
  //! \param geoBlk Running geometry block counter
  //! \param nBlock Running result block counter
  virtual bool writeGlvT(VTF* vtf, int iStep, int& geoBlk, int& nBlock) const;
  //! \brief Returns whether there are any pressure values to write to VTF.
  virtual bool hasTractionValues() const { return !presVal.empty(); }

protected:
  // Finite element quantities, i.e., indices into element matrices and vectors.
  // These indices will be identical for all elements in a model and can thus
  // be stored here, even when doing multi-threading. Note that the indices are
  // 1-based, since the value zero is used to signal non-existing matrix/vector.
  unsigned short int eK; //!< Index to element stiffness matrix
  unsigned short int eM; //!< Index to element mass matrix
  unsigned short int eS; //!< Index to element load vector

  // Physical properties
  Material* material;  //!< Material data and constitutive relation
  double    thickness; //!< Plate/shell thickness
  double    gravity;   //!< Gravitation constant

  LocalSystem* locSys;  //!< Local coordinate system for result output
  RealFunc*    presFld; //!< Pointer to pressure field

  mutable std::vector<Vec3Pair> presVal; //!< Pressure field point values
};

#endif
