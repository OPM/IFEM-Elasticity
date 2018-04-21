// $Id$
//==============================================================================
//!
//! \file main.C
//!
//! \date Apr 21 2018
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Main program for the isogeometric shell solver.
//!
//==============================================================================

#include "IFEM.h"
#include "SIMShell.h"
#include "NonlinearDriver.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "Utilities.h"
#include "Profiler.h"
#include <stdlib.h>
#include <string.h>


/*!
  \brief Main program for the isogeometric shell solver.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -hdf5 : Write primary and projected secondary solution to HDF5 file
  \arg -saveInc \a dtSave : Time increment between each result save to VTF/HDF5
  \arg -outPrec \a nDigit : Number of digits in solution component printout
  \arg -ztol \a eps : Zero tolerance for printing of solution norms
  \arg -stopTime \a t : Run simulation only up to specified stop time
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  int outPrec = 3;
  bool check = false;
  double zero_tol = 1.0e-8;
  double stopTime = 0.0;
  char* infile = nullptr;

  IFEM::Init(argc,argv,"Shell solver");

  for (int i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-outPrec") && i < argc-1)
      outPrec = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-ztol") && i < argc-1)
      zero_tol = atof(argv[++i]);
    else if (!strcmp(argv[i],"-stopTime") && i < argc-1)
      stopTime = atof(argv[++i]);
    else if (!strcmp(argv[i],"-check"))
      check = true;
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-nGauss <n>] [-check] [-hdf5]\n"
              <<"       [-vtf <format> [-nviz <nviz>] [-nu <nu>] [-nv <nv>]]\n"
              <<"       [-saveInc <dtSave>] [-outPrec <nd>] [-stopTime <t>]\n";
    return 0;
  }

  IFEM::getOptions().discretization = ASM::SplineC1;
  IFEM::cout <<"\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
  if (outPrec > 3)
    IFEM::cout <<"\nNorm- and component output precision: "<< outPrec;
  if (zero_tol != 1.0e-8)
    IFEM::cout <<"\nNorm output zero tolerance: "<< zero_tol;
  if (stopTime > 0.0)
    IFEM::cout <<"\nSimulation stop time: "<< stopTime;
  IFEM::cout << std::endl;

  utl::profiler->stop("Initialization");
  utl::profiler->start("Model input");

  // Create the shell FE model
  SIMShell model;
  NonlinearDriver simulator(model);

  // Read in solver and model definitions
  if (!simulator.read(infile))
    return 1;

  // Let the stop time specified on command-line override input file setting
  if (stopTime > 0.0)
    simulator.setStopTime(stopTime);

  model.opt.print(IFEM::cout,true) << std::endl;

  utl::profiler->stop("Model input");

  // Preprocess the model and establish data structures for the algebraic system
  if (!model.preprocess())
    return 2;

  if (model.opt.format >= 0)
  {
    // Save FE model to VTF file for visualization
    model.opt.nViz[2] = 1;
    if (!simulator.saveModel(infile))
      return 3;
  }

  if (check) return 0; // Data check only, no simulation

  // Initialize the linear equation solver and solution vectors
  simulator.initSol();

  DataExporter* writer = nullptr;
  if (model.opt.dumpHDF5(infile))
  {
    // Open HDF5 result database
    const std::string& fileName = model.opt.hdf5;
    IFEM::cout <<"\nWriting HDF5 file "<< fileName <<".hdf5"<< std::endl;

    writer = new DataExporter(true);
    writer->registerField("u","solution",DataExporter::SIM,
                          DataExporter::PRIMARY);
    writer->setFieldValue("u",&model,&simulator.getSolution());
    writer->registerWriter(new HDF5Writer(fileName,model.getProcessAdm()));
    writer->registerWriter(new XMLWriter(fileName,model.getProcessAdm()));
  }

  // Now invoke the main solution driver
  int status = simulator.solveProblem(writer,nullptr,0.0,zero_tol,outPrec);

  delete writer;
  return status;
}
