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
#include "Profiler.h"
#include "tinyxml2.h"
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
  virtual bool parse(const tinyxml2::XMLElement* elem)
  {
    if (!strcasecmp(elem->Value(),"newmarksolver"))
    {
      utl::getAttribute(elem,"initacc",doInitAcc);
      const tinyxml2::XMLElement* child = elem->FirstChildElement();
      for (; child; child = child->NextSiblingElement())
        params.parse(child);
    }
    else if (!strcasecmp(elem->Value(),"postprocessing"))
    {
      const tinyxml2::XMLElement* res = elem->FirstChildElement("resultpoints");
      if (res)
        // The point file is used only for point and line output (not for grid)
        if (res->FirstChildElement("point") || res->FirstChildElement("line"))
          utl::getAttribute(res,"file",rptFile);
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
    bool doProject = pi != Newmark::opt.project.end();
    if (doProject) proSol.resize(1);

    double nextSave = params.time.t + Newmark::opt.dtSave;

    // Open output file for result point print, if requested
    std::ostream* os = !rptFile.empty() ? new std::ofstream(rptFile) : nullptr;
    utl::LogStream* log = os ? new utl::LogStream(*os) : &IFEM::cout;
    std::streamsize rptPrec = outPrec > 0 ? outPrec : 3;

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
        Matrix projs(proSol.front());
        if (!Newmark::model.project(projs,this->realSolution(),
                                    pi->first,params.time))
          status += 6;
      }

      utl::profiler->start("Postprocessing");

      // Print solution components at the user-defined points
      Newmark::model.setMode(SIM::RECOVERY);
      this->dumpResults(params.time.t,*log,rptPrec,!os);

      if (params.hasReached(nextSave))
      {
        ++iStep;

        // Save solution variables to VTF
        if (Newmark::opt.format >= 0)
          status += this->saveStep(iStep,pi->second);

        // Save solution variables to HDF5
        if (writer && !writer->dumpTimeLevel(&params))
          status += 15;

        // Save solution variables to grid files, if specified
        if (!Newmark::model.saveResults(this->realSolutions(true),
                                        params.time.t,iStep))
          status += 16;

        nextSave = params.time.t + Newmark::opt.dtSave;
        if (nextSave > params.stopTime)
          nextSave = params.stopTime; // Always save the final step
      }

      // Save solution state to restart file
      if (restart && restart->dumpStep(params))
      {
        HDF5Restart::SerializeData data;
        if (this->serialize(data) && !restart->writeData(data))
          status += 17;
      }

      utl::profiler->stop("Postprocessing");
    }

    if (os)
    {
      delete log;
      delete os;
    }

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

  //! \brief Returns a pointer to the projected solution.
  const Vector* getProjection() const
  {
    return proSol.empty() ? nullptr : proSol.data();
  }

  //! \brief Dummy method required for template instantiation.
  Matrix* getNorms() { return nullptr; }

  //! \brief Overrides the stop time that was read from the input file.
  void setStopTime(double t) { params.stopTime = t; }

protected:
  using Newmark::saveStep;
  //! \brief Saves the converged solution to VTF file.
  int saveStep(int iStep, const std::string& prefix)
  {
    if (!this->saveStep(iStep,params.time.t))
      return 11;

    int ip = this->numSolution() - 2;
    if (ip > 0)
    {
      int nf = Newmark::model.getNoFields();
      if (!Newmark::model.writeGlvS1(this->realSolution(ip),iStep,
                                     Newmark::nBlock,params.time.t,
                                     "velocity",40,10+nf))
        return 12;

      if (!Newmark::model.writeGlvS1(this->realSolution(++ip),iStep,
                                     Newmark::nBlock,params.time.t,
                                     "acceleration",50,10+nf))
        return 13;
    }

    if (!proSol.empty())
      if (!Newmark::model.writeGlvP(proSol.front(),iStep,
                                    Newmark::nBlock,110,
                                    prefix.c_str()))
        return 14;

    return 0;
  }

private:
  std::string rptFile;   //!< Name of output file for point/line results
  bool        doInitAcc; //!< If \e true, calculate initial accelerations

  Vectors proSol; //!< Projected secondary solution

protected:
  TimeStep params; //!< Time stepping parameters
};

#endif
