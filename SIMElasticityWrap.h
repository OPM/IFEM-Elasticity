// $Id$
//==============================================================================
//!
//! \file SIMElasticityWrap.h
//!
//! \date April 21 2016
//!
//! \author Knut Morten Okstad
//!
//! \brief Wrapper equipping the elasticity solver with an ISolver interface.
//!
//==============================================================================

#ifndef _SIM_ELASTICITY_WRAP_H_
#define _SIM_ELASTICITY_WRAP_H_

#include "SIMElasticity.h"
#include "SIMsolution.h"


class DataExporter;


/*!
  \brief Driver wrapping an elasticity solver with an ISolver interface.
  \details The purpose of this class is to extend the SIMElasticity class with
  the required methods such that it fits the ISolver interface and thus can be
  used as class template argument to a SIMSolver instance for coupled solvers.
*/

template<class Dim>
class SIMElasticityWrap : public SIMElasticity<Dim>, public SIMsolution
{
protected:
  //! \brief The default constructor is protected as this is an interface class.
  SIMElasticityWrap();

public:
  //! \brief Empty destructor.
  virtual ~SIMElasticityWrap() {}

  //! \brief Registers solution fields for data output.
  //! \param exporter Result export handler
  void registerFields(DataExporter& exporter);

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param geoBlk Running geometry block counter
  bool saveModel(char* fileName, int& geoBlk, int&);

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param nBlock Running result block counter
  virtual bool saveStep(const TimeStep& tp, int& nBlock);

  //! \brief Serializes current internal state for restarting purposes.
  bool serialize(SerializeMap& data) const override;

  //! \brief Restores the internal state from serialized data.
  bool deSerialize(const SerializeMap& data) override;

  //! \brief Restores the basis from serialized data.
  bool deSerializeBasis(const SerializeMap& data);

  //! \brief Initializes the linear equation solver and solution vectors.
  //! \param[in] tp Time stepping parameters
  //! \param[in] withRF If \e true, reaction forces will be calculated
  virtual bool init(const TimeStep& tp, bool withRF = false);

  //! \brief Advances the time step one step forward.
  bool advanceStep(TimeStep& tp) override;

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp) = 0;

  // Due to the multiple inheritance, the compiler needs to be told which
  // version of this method to use (even though they have different signature)
  using SIMsolution::getSolution;

  //! \brief Overrides the parent class method to do nothing.
  bool postProcessNorms(Vectors&, Matrix*) override { return true; }
};

#endif
