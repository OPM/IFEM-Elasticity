// $Id$
//==============================================================================
//!
//! \file NLKirchhoffLoveShell.h
//!
//! \date Apr 20 2018
//!
//! \author Simen Skogholt Haave and Marit Gaarder Rakvaag / NTNU
//!
//! \brief Class for nonlinear Kirchhoff-Love thin shell problems.
//!
//==============================================================================

#ifndef _NL_KIRCHHOFF_LOVE_SHELL_H
#define _NL_KIRCHHOFF_LOVE_SHELL_H

#include "KirchhoffLoveShell.h"


/*!
  \brief Class representing the integrand of nonlinear thin shell problems.

  \details The formulation is based on Kirchhoff-Love shell theory
  and therefore requires second-derivatives of the basis functions.
*/

class NLKirchhoffLoveShell : public KirchhoffLoveShell
{
public:
  //! \brief Default constructor.
  NLKirchhoffLoveShell() {}
  //! \brief Empty destructor,
  virtual ~NLKirchhoffLoveShell() {}

  //! \brief Defines the solution mode before the element assembly is started.
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const;

  using KirchhoffLoveShell::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

private:
  //! \brief Evaluates the stiffness matrix and internal forces integrand.
  //! \param EK Element matrix to receive the stiffness contributions
  //! \param ES Element vector to receive the internal forc contributions
  //! \param[in] fe Finite element data at current point
  //! \param[in] G0 Co-variant basis vectors at the reference configuration
  //! \param[in] Gn Co-variant basis vectors at the actual configuration
  //! \param[in] H0 Hessian at the reference configuration
  //! \param[in] Hn Hessian at the actual configuration
  //! \param[in] X Cartesian coordinates of current point
  bool evalKandS(Matrix& EK, Vector& ES, const FiniteElement& fe,
                 const Matrix& G0, const Matrix& Gn,
                 const Matrix& H0, const Matrix& Hn, const Vec3& X) const;
};

#endif
