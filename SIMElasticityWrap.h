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
#include "SIMenums.h"
#include "SIMsolution.h"
#include "DataExporter.h"
#include "TimeStep.h"


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
  SIMElasticityWrap()
  {
    Dim::msgLevel = 1;
    Dim::myHeading = "Elasticity solver";
  }

public:
  //! \brief Empty destructor.
  virtual ~SIMElasticityWrap() {}

  //! \brief Registers solution fields for data output.
  //! \param exporter Result export handler
  void registerFields(DataExporter& exporter)
  {
    int flag = DataExporter::PRIMARY;
    if (!Dim::opt.pSolOnly)
      flag |= DataExporter::SECONDARY;
    exporter.registerField("u","solution",DataExporter::SIM,flag);
    exporter.setFieldValue("u",this,&this->getSolution());
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param geoBlk Running geometry block counter
  bool saveModel(char* fileName, int& geoBlk, int&)
  {
    if (Dim::opt.format < 0)
      return true;

    return this->writeGlvG(geoBlk,fileName);
  }

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param nBlock Running result block counter
  virtual bool saveStep(const TimeStep& tp, int& nBlock)
  {
    if (Dim::opt.format < 0 || tp.step%Dim::opt.saveInc > 0)
      return true;

    int iDump = 1 + tp.step/Dim::opt.saveInc;
    if (!this->writeGlvS(this->getSolution(),iDump,nBlock,tp.time.t,"u"))
      return false;

    return this->writeGlvStep(iDump,tp.time.t);
  }

  //! \brief Serializes current internal state for restarting purposes.
  virtual bool serialize(SerializeMap& data) const
  {
    if (!this->saveBasis(data) || !this->saveSolution(data,this->getName()))
      return false;

    data["Elasticity::Eext"] = SIMsolution::serialize(this->getExtEnerg(),1);

    return true;
  }

  //! \brief Restores the internal state from serialized data.
  virtual bool deSerialize(const SerializeMap& data)
  {
    if (!this->restoreSolution(data,this->getName()))
      return false;

    SerializeMap::const_iterator sit = data.find("Elasticity::Eext");
    if (sit != data.end())
      SIMsolution::deSerialize(sit->second,this->theExtEnerg(),1);

    return true;
  }

  //! \brief Restores the basis from serialized data.
  bool deSerializeBasis(const SerializeMap& data)
  {
    return this->restoreBasis(data);
  }

  //! \brief Initializes the linear equation solver and solution vectors.
  //! \param[in] tp Time stepping parameters
  //! \param[in] withRF If \e true, reaction forces will be calculated
  virtual bool init(const TimeStep& tp, bool withRF = false)
  {
    return (this->initSystem(Dim::opt.solver,1,1,0,withRF) &&
            this->initSolution(this->getNoDOFs(),this->getNoSolutions()) &&
            this->setMode(SIM::INIT) &&
            this->getIntegrand()->init(tp.time));
  }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp)
  {
    this->pushSolution(); // Update solution vectors between time steps
    return this->SIMElasticity<Dim>::advanceStep(tp);
  }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp) = 0;

  // Due to the multiple inheritance, the compiler needs to be told which
  // version of this method to use (even though they have different signature)
  using SIMsolution::getSolution;

  //! \brief Overrides the parent class method to do nothing.
  virtual bool postProcessNorms(Vectors&, Matrix*) { return true; }
};

#endif
