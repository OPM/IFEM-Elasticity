// $Id$
//==============================================================================
//!
//! \file ArcLengthDriver.h
//!
//! \date Apr 12 2018
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Arc-length solution driver for nonlinear isogeometric FEM simulators.
//!
//==============================================================================

#ifndef _ARC_LENGTH_DRIVER_H
#define _ARC_LENGTH_DRIVER_H

#include "NonlinearDriver.h"


/*!
  \brief Arc-length solution driver for nonlinear isogeometric FEM simulators.
  \details This class overrides the NonlinearDriver::solveStep() method,
  which does nonlinear quasi-static simulations with prescribed increment size,
  to implement a path-following solution procedure using the arc-length method.
  This can be used to solve nonlinear problems with instabilities, such as
  snap-through, snap-back, etc.
*/

class ArcLengthDriver : public NonlinearDriver
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  //! \param sim Reference to the spline FE model
  //! \param[in] adaptive If \e true, use adaptive mesh refinement
  explicit ArcLengthDriver(SIMbase& sim, bool = false, bool adaptive = false);

  //! \brief Advances the load step one step forward.
  virtual bool advanceStep(TimeStep& param, bool);

  //! \brief Solves the nonlinear equations by Newton-Raphson iterations.
  //! \param param Time stepping parameters
  //! \param[in] zero_tol Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual SIM::ConvStatus solveStep(TimeStep& param, SIM::SolutionMode,
                                    double zero_tol, std::streamsize outPrec);

  //! \brief Does nothing, overrides parent class method.
  virtual SIM::ConvStatus solveIteration(TimeStep&) { return SIM::FAILURE; }

protected:
  //! \brief Checks whether the nonlinear iterations have converged or diverged.
  virtual SIM::ConvStatus checkConvergence(TimeStep& param);
  //! \brief Does nothing, overrides parent class method.
  virtual bool lineSearch(TimeStep&) { return false; }

  using NonlinearDriver::parse;
  //! \brief Parses a data section from an XML document.
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] lambda Current load proportionality factor
  //! \param[in] iter Newton iteration counter
  //! \return Norm of the external load gradient, negative on failure
  double solveLinearizedSystem(double lambda, int iter = 0);

private:
  double arclen; //!< Arc-length parameter
  double beta;   //!< Arc-length parameter

  Vector tansol; //!< Tangential displacement vector
  Vector incsol; //!< Incremental solution of previous load step
};

#endif
