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
#include "DataExporter.h"


/*!
  \brief Driver wrapping an elasticity solver with an ISolver interface.
  \details The purpose of this class is to extend the SIMElasticity class with
  the required methods such that it fits the ISolver interface and thus can be
  used as class template argument to a SIMSolver instance for coupled solvers.
*/

template<class Dim> class SIMElasticityWrap : public SIMElasticity<Dim>
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

  //! \brief Initializes the solution vectors.
  //! \param[in] tp Time stepping parameters
  virtual bool init(const TimeStep& tp) = 0;

  //! \brief Computes the solution for the current time step.
  //! \param tp Time stepping parameters
  virtual bool solveStep(TimeStep& tp) = 0;

protected:
  //! \brief Returns a const reference to current solution vector.
  virtual const Vector& getSolution(int = 0) const = 0;
};

#endif
