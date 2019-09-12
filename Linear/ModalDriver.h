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


/*!
  \brief Driver for modal analysis of linear dynamic problems.
*/

class ModalDriver : public NewmarkDriver<NewmarkSIM>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit ModalDriver(SIMbase& sim) : NewmarkDriver<NewmarkSIM>(sim) {}
  //! \brief Empty destructor.
  virtual ~ModalDriver() {}

  //! \brief Calculates the current real solution vectors.
  virtual const Vectors& realSolutions();
  //! \brief Returns a const reference to the current real solution vector.
  virtual const Vector& realSolution(int i = 0) const;
  //! \brief Returns the number of solution vectors.
  virtual size_t numSolution() const;

  //! \brief Serialize solution state for restarting purposes.
  virtual bool serialize(HDF5Restart::SerializeData& data) const;
  //! \brief Set solution from a serialized state.
  virtual bool deSerialize(const HDF5Restart::SerializeData& data);

  //! \brief Dumps solution variables at user-defined points.
  virtual void dumpResults (double time, utl::LogStream& os,
                            std::streamsize precision, bool formatted) const;

  //! \brief Dumps the projected secondary solution for the eigenmodes.
  void dumpModes(utl::LogStream& os, std::streamsize precision) const;
};

#endif
