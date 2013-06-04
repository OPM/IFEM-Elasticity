// $Id$
//==============================================================================
//!
//! \file main_LinEl.C
//!
//! \date Jan 28 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Main program for the isogeometric linear elasticity solver.
//!
//==============================================================================

#include "IFEM.h"
#include "SIMLinEl.h"
#include "SIMLinElKL.h"
#include "SIMLinElBeamC1.h"
#include "ImmersedBoundaries.h"
#include "AdaptiveSIM.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "Utilities.h"
#include "VTF.h"
#include "Profiler.h"
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


/*!
  \brief Main program for the NURBS-based isogeometric linear elasticity solver.

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
  \arg -LR : Use LR-spline basis functions instead of tensorial splines/NURBS
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -hdf5 : Write primary and projected secondary solution to HDF5 file
  \arg -dumpASC : Dump model and solution to ASCII files for external processing
  \arg -ignore \a p1, \a p2, ... : Ignore these patches in the analysis
  \arg -eig \a iop : Eigenproblem solver to use (1...6)
  \arg -nev \a nev : Number of eigenvalues to compute
  \arg -ncv \a ncv : Number of Arnoldi vectors to use in the eigenvalue analysis
  \arg -shift \a shf : Shift value to use in the eigenproblem solver
  \arg -free : Ignore all boundary conditions (use in free vibration analysis)
  \arg -check : Data check only, read model and output to VTF (no solution)
  \arg -checkRHS : Check that the patches are modelled in a right-hand system
  \arg -vizRHS : Save the right-hand-side load vector on the VTF-file
  \arg -fixDup : Resolve co-located nodes by merging them into a single node
  \arg -2D : Use two-parametric simulation driver (plane stress)
  \arg -2Dpstrain : Use two-parametric simulation driver (plane strain)
  \arg -2Daxisymm : Use two-parametric simulation driver (axi-symmetric solid)
  \arg -KL : Use two-parametric simulation driver for Kirchhoff-Love plate
  \arg -Beam : Use one-parametric simulation driver for C1-continous beam
  \arg -adap : Use adaptive simulation driver with LR-splines discretization
  \arg -DGL2 : Estimate error using discrete global L2 projection
  \arg -CGL2 : Estimate error using continuous global L2 projection
  \arg -SCR : Estimate error using Superconvergent recovery at Greville points
  \arg -VDSA: Estimate error using Variational Diminishing Spline Approximations
  \arg -LSQ : Estimate error using through Least Square projections
  \arg -QUASI : Estimate error using Quasi-interpolation projections
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  SIMoptions dummy;
  std::vector<int> ignoredPatches;
  size_t adaptor = 0;
  int  i, iop = 0;
  bool checkRHS = false;
  bool vizRHS = false;
  bool fixDup = false;
  bool dumpASCII = false;
  bool twoD = false;
  bool KLp = false;
  bool Beam = false;
  char* infile = 0;

  int myPid = IFEM::Init(argc,argv);

  for (i = 1; i < argc; i++)
    if (dummy.parseOldOptions(argc,argv,i))
      ; // ignore the obsolete option
    else if (!strcmp(argv[i],"-dumpASC"))
      dumpASCII = myPid == 0; // not for parallel runs
    else if (!strcmp(argv[i],"-ignore"))
      while (i < argc-1 && isdigit(argv[i+1][0]))
	utl::parseIntegers(ignoredPatches,argv[++i]);
    else if (!strcmp(argv[i],"-vox") && i < argc-1)
      VTF::vecOffset[0] = atof(argv[++i]);
    else if (!strcmp(argv[i],"-voy") && i < argc-1)
      VTF::vecOffset[1] = atof(argv[++i]);
    else if (!strcmp(argv[i],"-voz") && i < argc-1)
      VTF::vecOffset[2] = atof(argv[++i]);
    else if (!strcmp(argv[i],"-plotSC"))
      SIMLinEl2D::GIpointsVTF = Immersed::plotCells = true;
    else if (!strcmp(argv[i],"-free"))
      SIMbase::ignoreDirichlet = true;
    else if (!strcmp(argv[i],"-check"))
      iop = 100;
    else if (!strcmp(argv[i],"-checkRHS"))
      checkRHS = true;
    else if (!strcmp(argv[i],"-vizRHS"))
      vizRHS = true;
    else if (!strcmp(argv[i],"-fixDup"))
      fixDup = true;
    else if (!strcmp(argv[i],"-Beam"))
      Beam = true;
    else if (!strcmp(argv[i],"-KL"))
      twoD = KLp = true;
    else if (!strncmp(argv[i],"-2Dpstra",8))
      twoD = SIMLinEl2D::planeStrain = true;
    else if (!strncmp(argv[i],"-2Daxi",6))
      twoD = SIMLinEl2D::axiSymmetry = true;
    else if (!strncmp(argv[i],"-2D",3))
      twoD = true;
    else if (!strncmp(argv[i],"-adap",5))
    {
      iop = 10;
      if (strlen(argv[i]) > 5)
        adaptor = atoi(argv[i]+5);
    }
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
	      <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n      "
	      <<" [-free] [-lag|-spec|-LR] [-2D[pstrain|axisymm]] [-KL|-Beam]"
	      <<" [-adap[<i>]] [-nGauss <n>]\n       [-vtf <format>"
	      <<" [-nviz <nviz>] [-nu <nu>] [-nv <nv>] [-nw <nw>]] [-hdf5]\n"
	      <<"       [-DGL2] [-CGL2] [-SCR] [-VDLSA] [-LSQ] [-QUASI]\n"
	      <<"       [-eig <iop> [-nev <nev>] [-ncv <ncv] [-shift <shf>]]\n"
	      <<"       [-ignore <p1> <p2> ...] [-fixDup]"
	      <<" [-checkRHS] [-check] [-dumpASC]\n";
    return 0;
  }

  if (myPid == 0)
  {
    const SIMoptions& opts = IFEM::getOptions();
    std::cout <<"\n >>> IFEM Linear Elasticity solver <<<"
	      <<"\n =====================================\n"
	      <<"\n Executing command:\n";
    for (i = 0; i < argc; i++) std::cout <<" "<< argv[i];
    std::cout <<"\n\nInput file: "<< infile
	      <<"\nEquation solver: "<< opts.solver
	      <<"\nNumber of Gauss points: "<< opts.nGauss[0];
    if (opts.format >= 0)
    {
      std::cout <<"\nVTF file format: "<< (opts.format ? "BINARY":"ASCII")
		<<"\nNumber of visualization points: "<< opts.nViz[0];
      if (!Beam)
      {
	std::cout <<" "<< opts.nViz[1];
	if (!twoD) std::cout <<" "<< opts.nViz[2];
      }
    }

    if (opts.eig > 0)
      std::cout <<"\nEigenproblem solver: "<< opts.eig
		<<"\nNumber of eigenvalues: "<< opts.nev
		<<"\nNumber of Arnoldi vectors: "<< opts.ncv
		<<"\nShift value: "<< opts.shift;
    if (opts.discretization == ASM::Lagrange)
      std::cout <<"\nLagrangian basis functions are used";
    else if (opts.discretization == ASM::Spectral)
      std::cout <<"\nSpectral basis functions are used";
    else if (opts.discretization == ASM::LRSpline)
      std::cout <<"\nLR-spline basis functions are used";
    if (SIMbase::ignoreDirichlet)
      std::cout <<"\nSpecified boundary conditions are ignored";
    if (fixDup)
      std::cout <<"\nCo-located nodes will be merged";
    if (checkRHS)
      std::cout <<"\nCheck that each patch has a right-hand coordinate system";
    if (!ignoredPatches.empty())
    {
      std::cout <<"\nIgnored patches:";
      for (size_t i = 0; i < ignoredPatches.size(); i++)
	std::cout <<" "<< ignoredPatches[i];
    }
    std::cout << std::endl;
  }
  utl::profiler->stop("Initialization");
  utl::profiler->start("Model input");

  // Create the simulation model
  SIMbase* model;
  if (Beam)
    model = new SIMLinElBeamC1();
  else if (KLp)
    model = new SIMLinElKL();
  else if (twoD)
    model = new SIMLinEl2D();
  else
    model = new SIMLinEl3D(checkRHS);

  SIMinput* theSim = model;
  AdaptiveSIM* aSim = 0;
  if (iop == 10)
  {
    theSim = aSim = new AdaptiveSIM(model);
    IFEM::getOptions().discretization = ASM::LRSpline;
  }
  else if (KLp && IFEM::getOptions().discretization != ASM::LRSpline)
    IFEM::getOptions().discretization = ASM::SplineC1;

  // Read in model definitions
  if (!theSim->read(infile))
    return 1;

  for (i = 1; i < argc; i++)
    if (!strcmp(argv[i],"-ignore"))
      while (i < argc-1 && isdigit(argv[i+1][0])) ++i;

  // Boundary conditions can be ignored only in generalized eigenvalue analysis
  if (model->opt.eig != 4 && model->opt.eig != 6)
    SIMbase::ignoreDirichlet = false;

  if (Beam) model->opt.nViz[1] = model->opt.nViz[2] = 1;
  if (twoD) model->opt.nViz[2] = 1;

  // Load vector visualization is not available when using additional viz-points
  if (model->opt.nViz[0] >2 || model->opt.nViz[1] >2 || model->opt.nViz[2] >2)
    vizRHS = false;

  if (myPid == 0)
  {
    std::cout <<"\n\nEquation solver: "<< model->opt.solver
	      <<"\nNumber of Gauss points: "<< model->opt.nGauss[0]
	      <<" "<< model->opt.nGauss[1];
    if (model->opt.format >= 0)
    {
      std::cout <<"\nVTF file format: "<< (model->opt.format ? "BINARY":"ASCII")
		<<"\nNumber of visualization points: "<< model->opt.nViz[0];
      if (!Beam)
      {
	std::cout <<" "<< model->opt.nViz[1];
	if (!twoD) std::cout <<" "<< model->opt.nViz[2];
      }
    }
    if (model->opt.eig > 0)
      std::cout <<"\nEigenproblem solver: "<< model->opt.eig
		<<"\nNumber of eigenvalues: "<< model->opt.nev
		<<"\nNumber of Arnoldi vectors: "<< model->opt.ncv
		<<"\nShift value: "<< model->opt.shift;
    if (model->opt.discretization == ASM::Lagrange)
      std::cout <<"\nLagrangian basis functions are used";
    else if (model->opt.discretization == ASM::Spectral)
      std::cout <<"\nSpectral basis functions are used";
    else if (model->opt.discretization == ASM::LRSpline)
      std::cout <<"\nLR-spline basis functions are used";
    std::cout << std::endl;

    model->printProblem(std::cout);
  }

  utl::profiler->stop("Model input");

  // Establish the FE data structures
  if (!model->preprocess(ignoredPatches,fixDup))
    return 1;

  SIMoptions::ProjectionMap& pOpt = model->opt.project;
  SIMoptions::ProjectionMap::const_iterator pit;

  // Set default projection method (tensor splines only)
  bool staticSol = iop + model->opt.eig%5 == 0 || iop == 10;
  if (model->opt.discretization < ASM::Spline || !staticSol)
    pOpt.clear(); // No projection if Lagrange/Spectral or no static solution
  else if (model->opt.discretization == ASM::Spline && pOpt.empty())
    pOpt[SIMoptions::GLOBAL] = "Greville point projection";

  if (model->opt.discretization < ASM::Spline && !model->opt.hdf5.empty())
  {
    std::cout <<"\n ** HDF5 output is available for spline discretization only."
	      <<" Deactivating...\n"<< std::endl;
    model->opt.hdf5.clear();
  }

  const char* prefix[pOpt.size()];
  if ((model->opt.format >= 0 || model->opt.dumpHDF5(infile)) && staticSol)
    for (i = 0, pit = pOpt.begin(); pit != pOpt.end(); pit++)
      prefix[i++] = pit->second.c_str();

  model->setQuadratureRule(model->opt.nGauss[0],true);

  Matrix eNorm, ssol;
  Vector displ, load;
  Vectors projs(pOpt.size()), gNorm;
  std::vector<Mode> modes;
  std::vector<Mode>::const_iterator it;
  int iStep = 1, nBlock = 0;

  if (aSim)
    aSim->setupProjections();

  DataExporter* exporter = NULL;
  if (model->opt.dumpHDF5(infile) && staticSol)
  {
    if (myPid == 0)
      std::cout <<"\nWriting HDF5 file "<< model->opt.hdf5
                <<".hdf5"<< std::endl;

    // Include secondary results only if no projection has been requested.
    // The secondary results will be projected anyway, but without the
    // nodal averaging across patch boundaries in case of multiple patches.
    int results = DataExporter::PRIMARY | 
                  DataExporter::NORMS   |
                  DataExporter::DISPLACEMENT;
    if (pOpt.empty()) results |= DataExporter::SECONDARY;

    exporter = new DataExporter(true);
    exporter->registerField("u","solution",DataExporter::SIM,results);
    exporter->setFieldValue("u",model, aSim ? &aSim->getSolution() : &displ);
    for (i = 0, pit = pOpt.begin(); pit != pOpt.end(); i++, pit++) {
      exporter->registerField(prefix[i], "projected", DataExporter::SIM,
                              DataExporter::SECONDARY, prefix[i]);
      exporter->setFieldValue(prefix[i], model,
                              aSim ? &aSim->getProjection(i) : &projs[i]);
    }
    exporter->registerWriter(new HDF5Writer(model->opt.hdf5));
    exporter->registerWriter(new XMLWriter(model->opt.hdf5));
    exporter->setNormPrefixes(prefix);
  }

  switch (iop+model->opt.eig) {
  case 0:
  case 5:
    // Static solution: Assemble [Km] and {R}
    model->setMode(SIM::STATIC);
    model->initSystem(model->opt.solver,1,1);
    model->setAssociatedRHS(0,0);
    if (!model->assembleSystem())
      return 2;
    else if (vizRHS)
      model->extractLoadVec(load);

    // Solve the linear system of equations
    if (!model->solveSystem(displ,1))
      return 3;

    // Project the FE stresses onto the splines basis
    model->setMode(SIM::RECOVERY);
    for (i = 0, pit = pOpt.begin(); pit != pOpt.end(); i++, pit++)
      if (!model->project(ssol,displ,pit->first))
	return 4;
      else
	projs[i] = ssol;

    if (myPid == 0 && !pOpt.empty())
      std::cout << std::endl;

    // Evaluate solution norms
    model->setQuadratureRule(model->opt.nGauss[1]);
    if (!model->solutionNorms(Vectors(1,displ),projs,eNorm,gNorm))
      return 4;

    if (myPid == 0)
    {
      std::cout <<"Energy norm |u^h| = a(u^h,u^h)^0.5   : "<< gNorm[0](1);
      std::cout	<<"\nExternal energy ((f,u^h)+(t,u^h)^0.5 : "<< gNorm[0](2);
      if (model->haveAnaSol() && gNorm[0].size() >= 4)
	std::cout <<"\nExact norm  |u|   = a(u,u)^0.5       : "<< gNorm[0](3)
		  <<"\nExact error a(e,e)^0.5, e=u-u^h      : "<< gNorm[0](4)
		  <<"\nExact relative error (%) : "
		  << gNorm[0](4)/gNorm[0](3)*100.0;
      size_t j = 1;
      for (pit = pOpt.begin(); pit != pOpt.end() && j < gNorm.size(); pit++,j++)
      {
	std::cout <<"\n\n>>> Error estimates based on "<< pit->second <<" <<<";
	std::cout <<"\nEnergy norm |u^r| = a(u^r,u^r)^0.5   : "<< gNorm[j](1);
	std::cout <<"\nError norm a(e,e)^0.5, e=u^r-u^h     : "<< gNorm[j](2);
	std::cout <<"\n- relative error (% of |u^r|) : "
		  << gNorm[j](2)/gNorm[j](1)*100.0;

        if (j == 0)
          continue;

	if (model->haveAnaSol())
	{
	  std::cout <<"\nExact error a(e,e)^0.5, e=u-u^r      : "<< gNorm[j](5)
		    <<"\n- relative error (% of |u|)   : "
		    << gNorm[j](5)/gNorm[0](3)*100.0;
	  std::cout <<"\nEffectivity index             : "
		    << gNorm[j](2)/gNorm[0](4);
	}

	std::cout <<"\nL2-norm |s^r| =(s^r,s^r)^0.5         : "<< gNorm[j](3);
	std::cout <<"\nL2-error (e,e)^0.5, e=s^r-s^h        : "<< gNorm[j](4);
        std::cout <<"\n- relative error (% of |s^r|) : "
                  << gNorm[j](4)/gNorm[j](3)*100.0;
      }
      std::cout << std::endl;

      model->dumpResults(displ,0.0,std::cout,true,6);
    }

    if (model->opt.eig == 0) break;

    // Linearized buckling: Assemble [Km] and [Kg]
    model->setMode(SIM::BUCKLING);
    model->initSystem(model->opt.solver,2,0);
    if (!model->assembleSystem(Vectors(1,displ)))
      return 5;

    // Solve the generalized eigenvalue problem
    if (!model->systemModes(modes))
      return 6;
    break;

  case 1:
  case 2:
    // Assemble and solve the regular eigenvalue problem
    model->setMode(SIM::STIFF_ONLY);
    model->initSystem(model->opt.solver,1,0);
    if (!model->assembleSystem())
      return 5;

    if (!model->systemModes(modes))
      return 6;
    break;

  case 10:
    // Adaptive simulation
    if (!aSim->initAdaptor(adaptor,4))
      break;

    if (exporter)
      exporter->setNormPrefixes(aSim->getNormPrefixes());

    while (true) {
      if (!aSim->solveStep(infile,iStep))
        return 5;
      else if (!aSim->writeGlv(infile,iStep,nBlock,4))
        return 6;
      else if (exporter)
        exporter->dumpTimeLevel(NULL,true);
      if (!aSim->adaptMesh(++iStep))
        break;
    }

  case 100:
    break; // Model check

  default:
    // Free vibration: Assemble [Km] and [M]
    model->setMode(SIM::VIBRATION);
    model->initSystem(model->opt.solver,2,0);
    if (!model->assembleSystem())
      return 5;

    if (!model->systemModes(modes))
      return 6;
  }

  utl::profiler->start("Postprocessing");

  if (iop != 10 && model->opt.format >= 0)
  {
    // Write VTF-file with model geometry
    if (!model->writeGlvG(nBlock,infile))
      return 7;

    // Write boundary tractions, if any
    if (!model->writeGlvT(iStep,nBlock))
      return 8;

    // Write Dirichlet boundary conditions
    if (!model->writeGlvBC(nBlock))
      return 8;

    // Write load vector to VTF-file
    if (!model->writeGlvV(load,"Load vector",iStep,nBlock))
      return 9;

    // Write solution fields to VTF-file
    if (!model->writeGlvS(displ,iStep,nBlock))
      return 10;

    // Write projected solution fields to VTF-file
    size_t i = 0;
    int iBlk = 100;
    for (pit = pOpt.begin(); pit != pOpt.end(); pit++, i++, iBlk += 10)
      if (!model->writeGlvP(projs[i],iStep,nBlock,iBlk,pit->second.c_str()))
	return 11;

    // Write eigenmodes
    bool isFreq = model->opt.eig==3 || model->opt.eig==4 || model->opt.eig==6;
    for (it = modes.begin(); it != modes.end(); it++)
      if (!model->writeGlvM(*it,isFreq,nBlock))
	return 12;

    // Write element norms
    if (!model->writeGlvN(eNorm,iStep,nBlock,prefix))
      return 13;

    model->writeGlvStep(1);
  }
  model->closeGlv();
  if (exporter && iop != 10)
    exporter->dumpTimeLevel();

  if (dumpASCII)
  {
    // Write (refined) model to g2-file
    std::ofstream osg(strcat(strtok(infile,"."),".g2"));
    osg.precision(18);
    std::cout <<"\nWriting updated g2-file "<< infile << std::endl;
    model->dumpGeometry(osg);
    if (!displ.empty())
    {
      // Write solution (control point values) to ASCII files
      std::ofstream osd(strcat(strtok(infile,"."),".dis"));
      osd.precision(18);
      std::cout <<"\nWriting deformation to file "<< infile << std::endl;
      model->dumpPrimSol(displ,osd,false);
      std::ofstream oss(strcat(strtok(infile,"."),".sol"));
      oss.precision(18);
      std::cout <<"\nWriting solution to file "<< infile << std::endl;
      model->dumpSolution(displ,oss);
    }
    if (!modes.empty())
    {
      // Write eigenvectors to ASCII files
      std::ofstream ose(strcat(strtok(infile,"."),".eig"));
      ose.precision(18);
      std::cout <<"\nWriting eigenvectors to file "<< infile << std::endl;
      for (it = modes.begin(); it != modes.end(); it++)
      {
	ose <<"# Eigenvector_"<< it->eigNo <<" Eigenvalue="<< it->eigVal <<"\n";
	model->dumpPrimSol(it->eigVec,ose,false);
      }
    }
  }

  utl::profiler->stop("Postprocessing");
  delete theSim;
  delete exporter;
  return 0;
}
