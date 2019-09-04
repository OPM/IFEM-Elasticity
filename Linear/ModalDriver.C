// $Id$
//==============================================================================
//!
//! \file ModalDriver.C
//!
//! \date Aug 29 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for modal analysis of linear dynamics problems.
//!
//==============================================================================

#include "NewmarkDriver.h"
#include "NewmarkSIM.h"


int modalSim (char* infile, size_t nM, SIMoutput* model, DataExporter* exporter)
{
  NewmarkDriver<NewmarkSIM> simulator(*model);

  // Read time integration setup
  if (!simulator.read(infile))
    return 1;

  // Initialize the modal time-domain simulation
  simulator.initPrm();
  simulator.initSolution(nM,3);
  simulator.printProblem();

  // Initialize the linear equation system.
  // Actually, we don't need any system matrices here
  // since we are only integrating the external load vector in space.
  if (!model->initSystem(LinAlg::DENSE,0,1))
    return 3;

  // Run the modal time integration
  return simulator.solveProblem(exporter,nullptr);
}
