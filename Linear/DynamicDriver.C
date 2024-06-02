// $Id$
//==============================================================================
//!
//! \file DynamicDriver.C
//!
//! \date May 31 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for linear dynamics problems.
//!
//==============================================================================

#include "IFEM.h"
#include "SIMoutput.h"
#include "NewmarkSIM.h"
#include "NewmarkDriver.h"
#include "ElasticityUtils.h"
#include "DataExporter.h"
#include "HDF5Writer.h"
#include "Profiler.h"


int dynamicSim (char* infile, SIMoutput* model, bool fixDup,
                double zero_tol, std::streamsize outPrec)
{
  IFEM::cout <<"\nUsing the linear dynamics simulation driver."<< std::endl;
  NewmarkDriver<NewmarkSIM> simulator(*model);

  // Read in solver and model definitions
  if (!simulator.read(infile))
    return 1;

  // Let the stop time specified on command-line override input file setting
  if (Elastic::time > 1.0)
    simulator.setStopTime(Elastic::time);

  model->opt.print(IFEM::cout,true) << std::endl;
  simulator.printProblem();

  utl::profiler->stop("Model input");

  // Preprocess the model and establish data structures for the algebraic system
  if (!model->preprocess({},fixDup))
    return 2;

  // Save FE model to VTF file for visualization
  if (model->opt.format >= 0 && !simulator.saveModel(infile))
    return 3;

  // Define the initial configuration and initialize the solution vectors
  simulator.initPrm();
  simulator.initSol(3);

  // Initialize the linear equation solver
  if (!simulator.initEqSystem(false,0))
    return 3;

  // Open HDF5 result database, if requested
  DataExporter* writer = nullptr;
  if (model->opt.dumpHDF5(infile))
  {
    const std::string& fileName = model->opt.hdf5;
    IFEM::cout <<"\nWriting HDF5 file "<< fileName <<".hdf5"<< std::endl;

    writer = new DataExporter(true,model->opt.saveInc);
    writer->registerWriter(new HDF5Writer(fileName,model->getProcessAdm()));

    int results = DataExporter::PRIMARY | DataExporter::DISPLACEMENT;
    writer->registerField("u","solution",DataExporter::SIM,results);
    writer->setFieldValue("u",model,&simulator.getSolution());
    results = -DataExporter::PRIMARY;
    writer->registerField("v","velocity",DataExporter::SIM,results);
    writer->setFieldValue("v",model,&simulator.getVelocity());
    writer->registerField("a","acceleration",DataExporter::SIM,results);
    writer->setFieldValue("a",model,&simulator.getAcceleration());
  }

  // Now invoke the main solution driver
  int status = simulator.solveProblem(writer,nullptr,zero_tol,outPrec);
  delete writer;
  return status;
}
