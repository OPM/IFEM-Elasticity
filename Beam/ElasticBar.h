// $Id$
//==============================================================================
//!
//! \file ElasticBar.h
//!
//! \date Aug 10 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class representing a nonlinear elastic bar.
//!
//==============================================================================

#ifndef _ELASTIC_BAR_H
#define _ELASTIC_BAR_H

#include "ElasticBase.h"


/*!
  \brief Class representing the integrand of a nonlinear elastic bar.

  \details Nonlinear elastic bar in 2D or 3D space with constant axial strain.

  The following three different strain meassures can be used:<OL>
  <LI> \a strain_meassure = 'E' : Engineering strain
       \f[ \varepsilon = \frac{L-L_0}{L_0} \f] </LI>
  <LI> \a strain_meassure = 'G' : Green strain
       \f[ \varepsilon = \frac{L^2-L_0^2}{2L_0^2} \f] </LI>
  <LI> \a strain_meassure = 'L' : Logarithmic strain
       \f[ \varepsilon = \log\left(\frac{L}{L_0}\right) \f] </LI>
  </OL>

  Reference: Kjell Magne Mathisen: Lecture 12:
  Formulation of Geometrically Nonlinear FE, October 2014.
*/

class ElasticBar : public ElasticBase
{
public:
  //! \brief The constructor initializes some parameters.
  //! \param[in] strm Strain meassure option ('E', 'G' or 'L')
  //! \param[in] nd Number of primary unknowns per node
  //! \param[in] ns Number of solution vectors in core
  ElasticBar(char strm, unsigned short int nd = 2, unsigned short int ns = 1);
  //! \brief Empty destructor.
  virtual ~ElasticBar() {}

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Defines the bar stiffness.
  void setStiffness(double stiff) { stiffness = stiff; }
  //! \brief Defines the mass density of the bar (mass per unit length).
  void setMass(double mass) { lineMass = fabs(mass); lumpedMass = mass >= 0.0; }

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const
  { return NO_DERIVATIVES | ELEMENT_CORNERS; }

  using ElasticBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

  using ElasticBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3&) const;

private:
  //! \brief Evaluates the axial strain.
  //! \param[in] LoverL0 Current length over initial length ratio
  double getStrain(double LoverL0) const;

protected:
  char strain_meassure; //!< Option for which strain meassure to use

  double stiffness; //!< Bar stiffness
  double lineMass;  //!< Bar mass pr unit length
  bool lumpedMass;  //!< Option of using lumped mass or consistent mass
};

#endif
