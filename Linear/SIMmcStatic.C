// $Id$
//==============================================================================
//!
//! \file SIMmcStatic.C
//!
//! \date Feb 12 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Monolithic coupling of linear static simulators.
//!
//==============================================================================

#include "SIMmcStatic.h"
#include "SIMoutput.h"
#include "SIMenums.h"
#include "DataExporter.h"
#include "IFEM.h"


int SIMmcStatic::solveStatic (const char* inpfile,
                              DataExporter* exporter,
                              double zero_tol, int outPrec)
{
  Vector displ;
  if (exporter)
  {
    std::string fieldName("u0");
    for (SIMoutput* sim : mySims)
    {
      ++fieldName[1];
      exporter->setFieldValue(fieldName, sim, &displ);
    }
  }

  // Static solution: Assemble [Km] and {R}
  int substep = 20;
  size_t i, nSim = mySims.size();
  for (i = 0; i < nSim; i++)
  {
    mySims[i]->printHeading(substep);
    mySims[i]->setMode(SIM::STATIC);
    mySims[i]->setQuadratureRule(mySims[i]->opt.nGauss[0],true,true);
    if (i == 0)
      mySims[i]->initSystem(mySims[i]->opt.solver);
    else
      mySims[i]->initSystem(mySims.front());
    if (!mySims[i]->assembleSystem())
      return 4;
  }

  // Solve the global linear system of equations
  this->printHeading(substep);
  if (!mySims.front()->solveSystem(displ,1))
    return 5;

  // Print result point values, if any
  for (SIMoutput* sim : mySims)
    if (sim->hasResultPoints())
    {
      sim->printHeading(substep);
      double old_tol = utl::zero_print_tol;
      if (zero_tol > 0.0) utl::zero_print_tol = zero_tol;
      sim->setMode(SIM::RECOVERY);
      sim->dumpResults(displ,0.0,IFEM::cout,true,outPrec);
      utl::zero_print_tol = old_tol;
    }

  if (exporter)
    exporter->dumpTimeLevel();

  if (opt.format < 0)
    return 0; // No VTF output

  // Helper class for cleaning the VTF objects on termination
  class ClearVTF
  {
    std::vector<SIMoutput*>& ourSims;
  public:
    ClearVTF(std::vector<SIMoutput*>& sims) : ourSims(sims) {}
    ~ClearVTF() { for (SIMoutput* sim : ourSims) sim->setVTF(nullptr); }
  } vtfCleaner(mySims);

  // Write VTF-file with model geometry
  int geoBlk = 0, nBlock = 0;
  std::vector<int> startBlk(nSim,0);
  if (!mySims.front()->writeGlvG(geoBlk,inpfile))
    return 12;

  for (i = 1; i < nSim; i++)
  {
    startBlk[i] = geoBlk;
    mySims[i]->setVTF(mySims.front()->getVTF());
    if (!mySims[i]->writeGlvG(geoBlk,0.0,true))
      return 12;
  }

  // Write solution fields to VTF-file
  for (i = 0; i < nSim; i++)
  {
    mySims[i]->setMode(SIM::RECOVERY);
    mySims[i]->setStartGeo(startBlk[i]);
    if (!mySims[i]->writeGlvS(displ,1,nBlock))
      return 16;
  }

  mySims.front()->writeGlvStep(1,0.0,-1);
  mySims.front()->closeGlv();

  return 0;

}
