// $Id$
//==============================================================================
//!
//! \file ElasticCable.h
//!
//! \date Aug 10 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class representing a 3D nonlinear elastic cable.
//!
//==============================================================================

#ifndef _ELASTIC_CABLE_H
#define _ELASTIC_CABLE_H

#include "ElasticBar.h"


/*!
  \brief Class representing the integrand of a 3D nonlinear elastic cable.

  \details This class implements the cable element described in the following
  paper: S. B. Raknes, et. al.: "Isogeometric rotation-free bending-stabilized
  cables: Statics, dynamics, bending strips and coupling with shells", Comput.
  Methods Appl. Mech, Engrg., 263 (2013) 127-143. The implementation is based on
  Siv Bente Raknes' Matlab code for the same element.
*/

class ElasticCable : public ElasticBar
{
public:
  //! \brief Default constructor.
  //! \param[in] n Number of consequtive solution vectors to reside in core
  ElasticCable(unsigned short int n = 1) : ElasticBar('G',3,n),
    EA(stiffness), EI(0.0) {}
  //! \brief Empty destructor.
  virtual ~ElasticCable() {}

  //! \brief Defines the bending stiffness.
  void setBendingStiffness(double stiff) { EI = stiff; }

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const { return SECOND_DERIVATIVES; }

  using ElasticBar::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  using ElasticBar::evalSol;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const FiniteElement& fe, const Vec3& X,
                       const std::vector<int>& MNPC) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld = 2) const { return fld == 1 ? 3 : 2; }
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField2Name(size_t i, const char* prefix = nullptr) const;

private:
  double& EA; //!< Axial stiffness
  double  EI; //!< Bending stiffness
};

#endif
