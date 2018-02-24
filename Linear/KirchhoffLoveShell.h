// $Id$
//==============================================================================
//!
//! \file KirchhoffLoveShell.h
//!
//! \date Feb 25 2018
//!
//! \author ... and ... / NTNU
//!
//! \brief Class for linear Kirchhoff-Love thin shell problems.
//!
//==============================================================================

#ifndef _KIRCHHOFF_LOVE_SHELL_H
#define _KIRCHHOFF_LOVE_SHELL_H

#include "KirchhoffLove.h"


/*!
  \brief Class representing the integrand of thin shell problems.

  \details The formulation is based on Kirchhoff-Love shell theory
  and therefore requires second-derivatives of the basis functions.
*/

class KirchhoffLoveShell : public KirchhoffLove
{
public:
  //! \brief Default constructor.
  KirchhoffLoveShell() : KirchhoffLove(3) {}
  //! \brief Empty destructor,
  virtual ~KirchhoffLoveShell() {}

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  using KirchhoffLove::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  using KirchhoffLove::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  using KirchhoffLove::evalSol;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const FiniteElement& fe,
                       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] s The FE stress resultant values at current point
  //! \param[in] eV Element solution vector
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] toLocal If \e true, transform to local coordinates (if defined)
  bool evalSol(Vector& s, const Vector& eV,
               const FiniteElement& fe, const Vec3& X,
               bool toLocal = false) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld = 2) const { return fld < 2 ? 3 : 6; }
  //! \brief Returns the name of the primary solution field.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField1Name(size_t i, const char* prefix) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField2Name(size_t i, const char* prefix) const;
};

#endif
