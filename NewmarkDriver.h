// $Id$
//==============================================================================
//!
//! \file NewmarkDriver.h
//!
//! \date Jul 9 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Newmark driver for isogeometric analysis of elastodynamics problems.
//!
//==============================================================================

#ifndef _NEWMARK_DRIVER_H
#define _NEWMARK_DRIVER_H

#include "IFEM.h"
#include "SIMoutput.h"
#include "SIMenums.h"
#include "DataExporter.h"
#include "HDF5Restart.h"
#include "TimeStep.h"
#include "Utilities.h"
#include "tinyxml.h"
#include <fstream>


/*!
  \brief Driver for isogeometric FEM analysis of elastodynamic problems.
*/

template<class Newmark> class NewmarkDriver : public Newmark
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  //! \param sim Reference to the spline FE model
  explicit NewmarkDriver(SIMbase& sim) : Newmark(sim) { doInitAcc = false; }
  //! \brief Empty destructor.
  virtual ~NewmarkDriver() {}

protected:
  using Newmark::parse;
  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem)
  {
    if (!strcasecmp(elem->Value(),"newmarksolver"))
    {
      utl::getAttribute(elem,"initacc",doInitAcc);
      const TiXmlElement* child = elem->FirstChildElement();
      for (; child; child = child->NextSiblingElement())
        params.parse(child);
    }
    else if (!strcasecmp(elem->Value(),"postprocessing"))
    {
      const TiXmlElement* respts = elem->FirstChildElement("resultpoints");
      if (respts)
        utl::getAttribute(respts,"file",pointfile);
    }

    bool ok = this->Newmark::parse(elem);
    if (ok && !strcasecmp(elem->Value(),"newmarksolver"))
      IFEM::cout <<"\talpha1 = "<< Newmark::alpha1
                 <<"  alpha2 = "<< Newmark::alpha2
                 <<"\n\tbeta = "<< Newmark::beta
                 <<"  gamma = "<< Newmark::gamma << std::endl;

    return ok;
  }

public:
  //! \brief Invokes the main time stepping simulation loop.
  //! \param writer HDF5 results exporter
  //! \param restart HDF5 restart handler
  //! \param[in] ztol Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  int solveProblem(DataExporter* writer, HDF5Restart* restart,
                   utl::LogStream*, bool, double,
                   double ztol = 1.0e-8, std::streamsize outPrec = 0)
  {
    return this->solveProblem(writer,restart,ztol,outPrec);
  }

  //! \brief Invokes the main time stepping simulation loop.
  //! \param writer HDF5 results exporter
  //! \param restart HDF5 restart handler
  //! \param[in] ztol Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  int solveProblem(DataExporter* writer, HDF5Restart* restart,
                   double ztol = 1.0e-8, std::streamsize outPrec = 0)
  {
    // Calculate initial accelerations
    if (doInitAcc && !this->initAcc(ztol,outPrec))
      return 4;

    SIMoptions::ProjectionMap::const_iterator pi = Newmark::opt.project.begin();
    bool doProject  = pi != Newmark::opt.project.end();
    double nextSave = params.time.t + Newmark::opt.dtSave;

    std::streamsize ptPrec = outPrec > 0 ? outPrec : 3;
    std::ostream* os = nullptr;
    if (!pointfile.empty())
      os = new std::ofstream(pointfile.c_str());

    // Invoke the time-step loop
    int status = 0;
    for (int iStep = 0; status == 0 && this->advanceStep(params);)
    {
      // Solve the dynamic FE problem at this time step
      if (this->solveStep(params,SIM::DYNAMIC,ztol,outPrec) != SIM::CONVERGED)
      {
        status = 5;
        break;
      }

      if (doProject)
      {
        // Project the secondary results onto the spline basis
        Newmark::model.setMode(SIM::RECOVERY);
        if (!Newmark::model.project(proSol,Newmark::solution.front(),
                                    pi->first,params.time))
          status += 6;
      }

      // Print solution components at the user-defined points
      if (os)
      {
        utl::LogStream log(os);
        this->dumpResults(params.time.t,log,ptPrec,false);
      }
      else
        this->dumpResults(params.time.t,IFEM::cout,ptPrec,true);

      if (params.hasReached(nextSave))
      {
        // Save solution variables to VTF
        if (Newmark::opt.format >= 0)
          if (!this->saveStep(++iStep,params.time.t) ||
              !Newmark::model.writeGlvS1(this->getVelocity(),iStep,
                                         Newmark::nBlock,params.time.t,
                                         "velocity",20) ||
              !Newmark::model.writeGlvS1(this->getAcceleration(),iStep,
                                         Newmark::nBlock,params.time.t,
                                         "acceleration",30) ||
              !Newmark::model.writeGlvP(proSol,iStep,Newmark::nBlock,110,
                                        pi->second.c_str()))
            status += 7;

        // Save solution variables to HDF5
        if (writer && !writer->dumpTimeLevel(&params))
          status += 8;

        nextSave = params.time.t + Newmark::opt.dtSave;
        if (nextSave > params.stopTime)
          nextSave = params.stopTime; // Always save the final step
      }

      // Save solution state to restart file
      if (restart && restart->dumpStep(params))
      {
        HDF5Restart::SerializeData data;
        if (this->serialize(data) && !restart->writeData(params,data))
          status += 9;
      }
    }

    delete os;

    return status;
  }

  //! \brief Serialize solution state for restarting purposes.
  //! \param data Container for serialized data
  virtual bool serialize(HDF5Restart::SerializeData& data) const
  {
    return params.serialize(data) && this->Newmark::serialize(data);
  }

  //! \brief Set solution from a serialized state.
  //! \param[in] data Container for serialized data
  virtual bool deSerialize(const HDF5Restart::SerializeData& data)
  {
    return params.deSerialize(data) && this->Newmark::deSerialize(data);
  }

  //! \brief Accesses the projected solution.
  const Vector& getProjection() const { return proSol; }

  //! \brief Overrides the stop time that was read from the input file.
  void setStopTime(double t) { params.stopTime = t; }

private:
  TimeStep params; //!< Time stepping parameters
  Matrix   proSol; //!< Projected secondary solution

  std::string pointfile; //!< Name of output file for point results
  bool        doInitAcc; //!< If \e true, calculate initial accelerations
};

#endif
