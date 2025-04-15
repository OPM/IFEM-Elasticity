// $Id$
//==============================================================================
//!
//! \file ModalDriver.h
//!
//! \date Aug 29 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for modal analysis of linear dynamics problems.
//!
//==============================================================================

#ifndef _MODAL_DRIVER_H
#define _MODAL_DRIVER_H

#include "NewmarkDriver.h"
#include "NewmarkSIM.h"
#include "EigenModeSIM.h"


/*!
  \brief Driver for modal analysis of linear dynamic problems.
*/

class ModalDriver : public NewmarkDriver<NewmarkSIM>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit ModalDriver(SIMbase& sim, bool qs = false)
    : NewmarkDriver<NewmarkSIM>(sim) { qstatic = qs; }

  //! \brief Calculates/returns the current real solution vectors.
  virtual const Vectors& realSolutions(bool returnCurrent);
  //! \brief Returns a const reference to the current real solution vector.
  virtual const Vector& realSolution(int i = 0) const;
  //! \brief Returns the number of solution vectors.
  virtual size_t numSolution() const;

  //! \brief Serialize solution state for restarting purposes.
  virtual bool serialize(HDF5Restart::SerializeData& data) const;
  //! \brief Set solution from a serialized state.
  virtual bool deSerialize(const HDF5Restart::SerializeData& data);

  //! \brief Dumps solution variables at user-defined points.
  virtual void dumpResults(double time, utl::LogStream& os,
                           std::streamsize precision, bool formatted) const;

  //! \brief Dumps the projected secondary solution for the eigenmodes.
  void dumpModes(utl::LogStream& os, std::streamsize precision) const;

  //! \brief Checks whether the corrector iterations have converged or diverged.
  SIM::ConvStatus checkConvergence(TimeStep& tp);
  //! \brief Calculates predicted velocities and accelerations.
  virtual bool predictStep(TimeStep& tp);
  //! \brief Updates configuration variables (solution vector) in an iteration.
  virtual bool correctStep(TimeStep& tp, bool);

private:
  bool qstatic; //!< If \e true, use quasi-static simulation driver
};


/*!
  \brief Driver creating a time history from a set of eigenmode shapes.
*/

class ModesHistorySIM : public EigenModeSIM
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit ModesHistorySIM(SIMbase& sim) : EigenModeSIM(sim) {}

protected:
  using EigenModeSIM::parse;
  //! \brief Parses a data section from an XML document.
  virtual bool parse(const tinyxml2::XMLElement* elem);

public:
  //! \brief Initializes the solver and runs through the time history.
  int solve(char* infile, double ztol = 1.0e-8, std::streamsize outPrec = 0);

private:
  TimeStep params; //!< Time stepping parameters
};

#endif
