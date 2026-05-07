// $Id$
//==============================================================================
//!
//! \file main_NonLinEl.C
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Main program for the isogeometric finite deformation solver.
//!
//==============================================================================

#include "SIMFiniteDefEl.h"
#include "SIMLinEl.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "HHTSIM.h"
#include "GenAlphaSIM.h"
#include "NewmarkNLSIM.h"
#include "NewmarkDriver.h"
#include "ArcLengthDriver.h"
#include "HDF5Writer.h"
#include "Utilities.h"
#include "Profiler.h"
#include "VTF.h"
#include "NLargs.h"
#include <filesystem>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cctype>

#ifndef USE_OPENMP
extern std::vector<int> dbgElms; //!< List of elements for additional output
#endif


/*!
  \brief Reads the input file and invokes the main simulation driver.
*/

template<class Simulator>
int runSimulator (Simulator& simulator, SIMoutput* model, char* infile,
                  const std::vector<int>& ignoredPatches, bool fixDup,
                  char printMax, double dtDump, double stopTime,
                  double zero_tol, int outPrec, bool dumpNodeMap)
{
  utl::profiler->start("Model input");

  std::ostream* oss = nullptr;

  // Lambda function cleaning heap-allocated objects before exiting.
  auto&& exitSim = [oss,model](int status)
  {
    delete oss;
    delete model;
    return status;
  };

  // Read in solver and model definitions
  if (!simulator.read(infile))
    return exitSim(1);

  // Let the stop time specified on command-line override input file setting
  if (stopTime > 0.0)
    simulator.setStopTime(stopTime);

  model->opt.print(IFEM::cout,true) << std::endl;
  simulator.printProblem();

  utl::profiler->stop("Model input");

  // Preprocess the model and establish data structures for the algebraic system
  if (!model->preprocess(ignoredPatches,fixDup))
    return exitSim(2);

  if (model->opt.format >= 0)
  {
    // Save FE model to VTF file for visualization
    if (model->getNoSpaceDim() < 3)
      model->opt.nViz[2] = 1;
    if (!simulator.saveModel(infile))
      return exitSim(4);
    else if (stopTime < 0.0 && !model->writeGlvStep(1))
      return exitSim(4);
  }

  if (dtDump < 0.0)
  {
    // Write (refined?) model to g2-file
    strcat(strtok(infile,"."),".g2");
    IFEM::cout <<"\nWriting updated g2-file "<< infile << std::endl;
    std::ofstream osg(infile);
    model->dumpGeometry(osg);

    // Open ASCII file for solution dump
    if (stopTime >= 0.0)
    {
      strcat(strtok(infile,"."),".sol");
      oss = new std::ofstream(infile);
      *oss <<"#NPoints="<< model->getNoNodes() <<"\n";
    }
  }

  if (stopTime < 0.0)
    return exitSim(0); // model check

  size_t numPatch = 1;
  const Elasticity* lelp;
  if (!(lelp = dynamic_cast<const Elasticity*>(model->getProblem())))
    printMax = false;
  else if (printMax == 'P')
    numPatch = model->getFEModel().size();
  if (printMax)
    const_cast<Elasticity*>(lelp)->initMaxVals(numPatch);

  if (model->opt.discretization < ASM::Spline && !model->opt.hdf5.empty())
  {
    IFEM::cout <<"\n ** HDF5 output is available for spline discretization only"
               <<". Deactivating...\n"<< std::endl;
    model->opt.hdf5.clear();
  }

  // We only allow the version=2 global L2-projection here
  SIMoptions::ProjectionMap& pOpt = model->opt.project;
  if (pOpt.find(SIMoptions::CGL2) == pOpt.end() &&
      pOpt.find(SIMoptions::CGL2_INT) == pOpt.end())
    pOpt.clear();
  else
  {
    pOpt.clear();
    pOpt[SIMoptions::CGL2_INT] = "Global L2 projection";
  }
  SIMoptions::ProjectionMap::const_iterator pit = pOpt.begin();

  // Define the initial configuration
  NewmarkSIM* dynSim = dynamic_cast<NewmarkSIM*>(&simulator);
  simulator.initPrm();
  simulator.initSol(dynSim ? 3 : 2);

  // Initialize the linear equation solver
  if (!simulator.initEqSystem(!dynSim, dynSim ? 0 : model->getNoFields()))
    return exitSim(3);

  // Load solution state from serialized data in case of restart
  if (!simulator.checkForRestart())
    return exitSim(5);

  // Open HDF5 result database
  DataExporter* writer = nullptr;
  if (model->opt.dumpHDF5(infile))
  {
    const std::string& fileName = model->opt.hdf5;
    IFEM::cout <<"\nWriting HDF5 file "<< fileName <<".hdf5"<< std::endl;

    // Include secondary results only if no projection has been requested.
    // The secondary results will be projected anyway, but without the
    // nodal averaging across patch boundaries in case of multiple patches.
    int results = DataExporter::PRIMARY;
    if (pit == pOpt.end() && !model->opt.pSolOnly)
      results |= DataExporter::SECONDARY;
    if (dumpNodeMap)
      results |= DataExporter::L2G_NODE;
    if (model->opt.saveNorms)
      results |= DataExporter::NORMS;
    if (model->hasElementActivator())
      results |= DataExporter::ELEMENT_MASK;

    writer = new DataExporter(true,model->opt.saveInc);
    writer->registerWriter(new HDF5Writer(fileName,model->getProcessAdm()));
    writer->registerField("u","solution",DataExporter::SIM,results);
    writer->setFieldValue("u",model,&simulator.getSolution(),
                          nullptr,simulator.getNorms());
    if (dynSim)
    {
      writer->registerField("v","velocity",DataExporter::SIM,
                            -DataExporter::PRIMARY);
      writer->setFieldValue("v",model,&dynSim->getVelocity());
      writer->registerField("a","acceleration",DataExporter::SIM,
                            -DataExporter::PRIMARY);
      writer->setFieldValue("a",model,&dynSim->getAcceleration());
    }
    if (pit != pOpt.end())
    {
      writer->registerField("sigma","projected",DataExporter::SIM,
                            DataExporter::SECONDARY,pit->second.c_str());
      writer->setFieldValue("sigma",model,simulator.getProjection());
    }
  }

  HDF5Restart* restart = nullptr;
  if (model->opt.restartInc > 0)
  {
    std::string hdf5file(infile);
    if (!model->opt.hdf5.empty())
      hdf5file = model->opt.hdf5 + "_restart";
    else
      hdf5file.replace(hdf5file.find_last_of('.'),std::string::npos,"_restart");
    const size_t idot = hdf5file.size();
    for (int i = 1; std::filesystem::exists(hdf5file + ".hdf5"); i++)
      hdf5file = hdf5file.substr(0,idot) + std::to_string(i);
    IFEM::cout <<"\nWriting HDF5 file "<< hdf5file <<".hdf5"<< std::endl;
    restart = new HDF5Restart(hdf5file,model->getProcessAdm(),
                              model->opt.restartInc);
  }

  if (pit != pOpt.end())
    IFEM::cout <<"\n"<< pit->second <<" will be used to compute"
               <<"\nsmoothed secondary solution fields."<< std::endl;

  // Now invoke the main solution driver
  utl::LogStream log(oss);
  int status = simulator.solveProblem(writer, restart, oss ? &log : nullptr,
                                      printMax, std::abs(dtDump),
                                      zero_tol, outPrec);

  delete writer;
  delete restart;
  return exitSim(status);
}


/*!
  \brief Main program for the isogeometric finite deformation solver.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -lag : Use Lagrangian basis functions instead of splines/NURBS
  \arg -spec : Use Spectral basis functions instead of splines/NURBS
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -hdf5 : Write primary and projected secondary solution to HDF5 file
  \arg -printMax : Print out maximum point-wise stresses
  \arg -printMaxPatch : Print out patch-wise maximum point-wise stresses
  \arg -saveInc \a dtSave : Time increment between each result save to VTF/HDF5
  \arg -dumpInc \a dtDump [raw] : Time increment between each solution dump
  \arg -dumpNodMap : Dump Local-to-global node number mapping to HDF5
  \arg -outPrec \a nDigit : Number of digits in solution component printout
  \arg -ztol \a eps : Zero tolerance for printing of solution norms
  \arg -ignore \a p1, \a p2, ... : Ignore these patches in the analysis
  \arg -check : Data check only, read model and output to VTF (no solution)
  \arg -checkRHS : Check that the patches are modelled in a right-hand system
  \arg -fixDup : Resolve co-located nodes by merging them into a single node
  \arg -stopTime \a t : Run simulation only up to specified stop time
  \arg -2D : Use two-parametric simulation driver (plane stress)
  \arg -2Dpstrain : Use two-parametric simulation driver (plane strain)
  \arg -2Daxisymm : Use two-parametric simulation driver (axi-symmetric solid)
  \arg -UL : Use updated Lagrangian formulation with nonlinear material
  \arg -MX<pord> : Mixed formulation with internal discontinuous pressure
  \arg -mixed : Mixed formulation with continuous pressure and volumetric change
  \arg -Mixed : Same as -mixed, but use C^(p-1) continuous displacement basis
  \arg -Fbar<nvp> : Use the F-bar formulation
  \arg -linear : Do a linear analysis only (no iterations)
  \arg -free : Ignore all boundary conditions (use in dynamics analysis)
  \arg -adap : Use adaptive simulation driver with LR-splines discretization
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);

  std::vector<int> ignoredPatches;
  int outPrec = 3;
  double dtDump = 0.0;
  double zero_tol = 1.0e-8;
  double stopTime = 0.0;
  char* infile = nullptr;
  NLargs args;

  IFEM::Init(argc,argv,"Finite Deformation Nonlinear solver");

  for (int i = 1; i < argc; i++)
    if (argv[i] == infile || args.parseArg(argv[i]))
      ; // ignore the input file on the second pass
    else if (SIMoptions::ignoreOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-outPrec") && i < argc-1)
      outPrec = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-ztol") && i < argc-1)
      zero_tol = atof(argv[++i]);
    else if (!strcmp(argv[i],"-dumpInc") && i < argc-1)
    {
      dtDump = atof(argv[++i]);
      if (++i < argc && !strcmp(argv[i],"raw"))
	dtDump *= -1;
      else
	--i;
    }
#ifndef USE_OPENMP
    else if (!strcmp(argv[i],"-dbgElm"))
      while (i < argc-1 && isdigit(argv[i+1][0]))
	utl::parseIntegers(dbgElms,argv[++i]);
#endif
    else if (!strcmp(argv[i],"-ignore"))
      while (i < argc-1 && isdigit(argv[i+1][0]))
	utl::parseIntegers(ignoredPatches,argv[++i]);
    else if (!strcmp(argv[i],"-vox") && i < argc-1)
      VTF::vecOffset[0] = atof(argv[++i]);
    else if (!strcmp(argv[i],"-voy") && i < argc-1)
      VTF::vecOffset[1] = atof(argv[++i]);
    else if (!strcmp(argv[i],"-voz") && i < argc-1)
      VTF::vecOffset[2] = atof(argv[++i]);
    else if (!strcmp(argv[i],"-free"))
      SIMbase::ignoreDirichlet = true;
    else if (!strcmp(argv[i],"-stopTime") && i < argc-1)
      stopTime = atof(argv[++i]);
    else if (!strcmp(argv[i],"-check"))
      stopTime = -1.0;
    else if (!infile)
    {
      infile = argv[i];
      if (strcasestr(infile,".xinp"))
      {
        if (!args.readXML(infile,false))
          return 1;
        i = 0; // start over and let command-line options override input file
      }
    }
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
	      <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
	      <<"       [-lag|-spec] [-2D[pstrain|axis]] [-nGauss <n>]\n"
	      <<"       [-UL|-MX[<p>]|-[M|m]ixed|-Fbar<nvp>]\n"
	      <<"       [-linear] [-adap] [-arclen|-HHT|-GA] [-free]\n"
	      <<"       [-hdf5 [<filename>] [-dumpNodeMap]]\n"
	      <<"       [-vtf <format> [-nviz <nviz>]"
	      <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]]\n      "
	      <<" [-saveInc <dtSave>] [-dumpInc <dtDump> [raw]]"
	      <<" [-outPrec <nd>]\n       [-ztol <eps>] [-ignore <p1> <p2> ...]"
	      <<" [-fixDup] [-checkRHS] [-check]\n"
	      <<"       [-printMax[Patch]] [-stopTime <t>]\n";
    return 0;
  }

  if (IFEM::getOptions().discretization == ASM::Spline && args.adaptive)
    IFEM::getOptions().discretization = ASM::LRSpline;
  else if (IFEM::getOptions().discretization < ASM::LRSpline)
    args.adaptive = false;

  IFEM::cout <<"\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
  if (dtDump > 0.0)
    IFEM::cout <<"\nTime between each primary solution dump: "<< dtDump;
  if (SIMbase::ignoreDirichlet)
    IFEM::cout <<"\nSpecified boundary conditions are ignored";
  if (args.fixDup)
    IFEM::cout <<"\nCo-located nodes will be merged";
  if (args.checkRHS)
    IFEM::cout <<"\nCheck that each patch has a right-hand coordinate system";
  if (!ignoredPatches.empty())
  {
    IFEM::cout <<"\nIgnored patches:";
    for (int ip : ignoredPatches) IFEM::cout <<" "<< ip;
  }
  if (outPrec != 3)
    IFEM::cout <<"\nNorm- and component output precision: "<< outPrec;
  else if (args.algor > STATIC)
    outPrec = 0;
  if (zero_tol != 1.0e-8)
    IFEM::cout <<"\nNorm output zero tolerance: "<< zero_tol;
  IFEM::cout << std::endl;

  utl::profiler->start("Model input");

  bool linear = !args.options.empty() && args.options.front() == SIM::LINEAR;

  SIMoutput* model;
  if (linear && args.algor > STATIC)
  {
    // Create the linear continuum model
    if (args.twoD)
      model = new SIMLinEl<SIM2D>(args.checkRHS);
    else
      model = new SIMLinEl<SIM3D>(args.checkRHS);

    if (args.algor == GENALPHA)
    {
      // Invoke the linear generalized alpha time integration
      NewmarkDriver<GenAlphaSIM> simulator(*model);
      return runSimulator(simulator,model,infile,ignoredPatches,args.fixDup,
                          args.printMax,dtDump,stopTime,zero_tol,outPrec,
                          args.dNodeMap);
    }

    // Invoke the linear Newmark time integration
    NewmarkDriver<NewmarkSIM> simulator(*model);
    return runSimulator(simulator,model,infile,ignoredPatches,args.fixDup,
                        args.printMax,dtDump,stopTime,zero_tol,outPrec,
                        args.dNodeMap);
  }

  // Create the nonlinear continuum model
  if (args.twoD)
    model = new SIMFiniteDefEl<SIM2D>(args.checkRHS,args.options);
  else
    model = new SIMFiniteDefEl<SIM3D>(args.checkRHS,args.options);

  switch (args.algor) {
  case STATIC:
  {
    // Invoke the nonlinear quasi-static solver with fixed load increments
    NonlinearDriver simulator(*model,linear,args.adaptive);
    return runSimulator(simulator,model,infile,ignoredPatches,args.fixDup,
                        args.printMax,dtDump,stopTime,zero_tol,outPrec,
                        args.dNodeMap);
  }
  case ARCLEN:
  {
    // Invoke the nonlinear quasi-static arc-length solver
    ArcLengthDriver simulator(*model,args.adaptive);
    return runSimulator(simulator,model,infile,ignoredPatches,args.fixDup,
                        args.printMax,dtDump,stopTime,zero_tol,outPrec,
                        args.dNodeMap);
  }
  case NEWHHT:
  {
    // Invoke the nonlinear HHT time integration
    NewmarkDriver<HHTSIM> simulator(*model);
    return runSimulator(simulator,model,infile,ignoredPatches,args.fixDup,
                        args.printMax,dtDump,stopTime,zero_tol,outPrec,
                        args.dNodeMap);
  }
  case OLDHHT:
  case GENALPHA:
  {
    // Invoke the nonlinear Newmark time integration
    NewmarkDriver<NewmarkNLSIM> simulator(*model);
    return runSimulator(simulator,model,infile,ignoredPatches,args.fixDup,
                        args.printMax,dtDump,stopTime,zero_tol,outPrec,
                        args.dNodeMap);
  }
  default:
    return -1; // Unknown driver
  }
}
