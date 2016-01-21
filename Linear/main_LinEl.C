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
#include "SIMElasticBar.h"
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
  \arg -2DKL : Use two-parametric simulation driver for Kirchhoff-Love plate
  \arg -1D : Use one-parametric simulation driver for beam with rotational DOFs
  \arg -1DKL : Use one-parametric simulation driver for C1-continous beam
  \arg -1DC1 : Use one-parametric simulation driver for C1-continous cable
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

  std::vector<int> ignoredPatches;
  size_t adaptor = 0;
  int  i, iop = 0;
  bool checkRHS = false;
  bool vizRHS = false;
  bool fixDup = false;
  bool dumpASCII = false;
  bool twoD = false;
  bool KLp = false;
  bool oneD = false;
  bool isC1 = false;
  bool noProj = false;
  bool noError = false;
  char* infile = NULL;
  Elasticity::wantPrincipalStress = true;

  int myPid = IFEM::Init(argc,argv);

  for (i = 1; i < argc; i++)
    if (SIMoptions::ignoreOldOptions(argc,argv,i))
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
    else if (!strcmp(argv[i],"-1DC1"))
      oneD = isC1 = true;
    else if (!strcmp(argv[i],"-1DKL"))
      oneD = isC1 = KLp = true;
    else if (!strncmp(argv[i],"-1D",3))
      oneD = true;
    else if (!strcmp(argv[i],"-2DKL"))
      twoD = isC1 = KLp = true;
    else if (!strncmp(argv[i],"-2Dpstra",8))
      twoD = SIMLinEl2D::planeStrain = true;
    else if (!strncmp(argv[i],"-2Daxi",6))
      twoD = SIMLinEl2D::axiSymmetry = true;
    else if (!strncmp(argv[i],"-2D",3))
      twoD = true;
    else if (!strncmp(argv[i],"-noP",4))
      noProj = true;
    else if (!strncmp(argv[i],"-noE",4))
      noError = true;
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
              <<" <inputfile> [-dense|-spr|-superlu[<nt>]|-samg|-petsc]\n"
              <<"       [-lag|-spec|-LR] [-1D[C1|KL]|-2D[pstrain|axisymm|KL]]"
              <<" [-nGauss <n>]\n       [-hdf5] [-vtf <format> [-nviz <nviz>]"
              <<" [-nu <nu>] [-nv <nv>] [-nw <nw>]]\n       [-adap[<i>]]"
              <<" [-DGL2] [-CGL2] [-SCR] [-VDLSA] [-LSQ] [-QUASI]\n      "
              <<" [-eig <iop> [-nev <nev>] [-ncv <ncv] [-shift <shf>] [-free]]"
              <<"\n       [-ignore <p1> <p2> ...] [-fixDup]"
              <<" [-checkRHS] [-check] [-dumpASC]\n";
    return 0;
  }

  if (iop == 10)
    IFEM::getOptions().discretization = ASM::LRSpline;
  else if (isC1 && IFEM::getOptions().discretization != ASM::LRSpline)
    IFEM::getOptions().discretization = ASM::SplineC1;

  IFEM::cout <<"\n >>> IFEM Linear Elasticity solver <<<"
             <<"\n =====================================\n"
             <<"\n Executing command:\n";
  for (i = 0; i < argc; i++) IFEM::cout <<" "<< argv[i];
  IFEM::cout <<"\n\nInput file: "<< infile;
  IFEM::getOptions().print(IFEM::cout);
  if (SIMbase::ignoreDirichlet)
    IFEM::cout <<"\nSpecified boundary conditions are ignored";
  if (fixDup)
    IFEM::cout <<"\nCo-located nodes will be merged";
  if (checkRHS && !oneD && !KLp)
    IFEM::cout <<"\nCheck that each patch has a right-hand coordinate system";
  if (!ignoredPatches.empty())
  {
    IFEM::cout <<"\nIgnored patches:";
    for (size_t i = 0; i < ignoredPatches.size(); i++)
      IFEM::cout <<" "<< ignoredPatches[i];
  }
  IFEM::cout << std::endl;

  utl::profiler->stop("Initialization");
  utl::profiler->start("Model input");

  // Create the simulation model
  SIMoutput* model;
  if (oneD)
  {
    if (KLp)
      model = new SIMLinElBeamC1();
    else if (IFEM::getOptions().eig == 4 || IFEM::getOptions().eig == 6)
      model = new SIMElasticDynBar();
    else
      model = new SIMElasticBar();
  }
  else if (KLp)
    model = new SIMLinElKL();
  else if (twoD)
    model = new SIMLinEl2D(checkRHS);
  else
    model = new SIMLinEl3D(checkRHS);

  SIMinput* theSim = model;
  AdaptiveSIM* aSim = NULL;
  if (iop == 10)
    theSim = aSim = new AdaptiveSIM(model);

  // Read in model definitions
  if (!theSim->read(infile))
    return 1;

  // Boundary conditions can be ignored only in generalized eigenvalue analysis
  if (model->opt.eig != 4 && model->opt.eig != 6)
    SIMbase::ignoreDirichlet = false;

  if (oneD) model->opt.nViz[1] = model->opt.nViz[2] = 1;
  if (twoD) model->opt.nViz[2] = 1;

  // Load vector visualization is not available when using additional viz-points
  if (model->opt.nViz[0] >2 || model->opt.nViz[1] >2 || model->opt.nViz[2] >2)
    vizRHS = false;

  model->opt.print(IFEM::cout,true) << std::endl;

  utl::profiler->stop("Model input");

  // Establish the FE data structures
  if (!model->preprocess(ignoredPatches,fixDup))
    return 1;

  SIMoptions::ProjectionMap& pOpt = model->opt.project;
  SIMoptions::ProjectionMap::const_iterator pit;

  // Set default projection method (tensor splines only)
  bool staticSol = iop + model->opt.eig%5 == 0 || iop == 10;
  if (model->opt.discretization < ASM::Spline || !staticSol || noProj)
    pOpt.clear(); // No projection if Lagrange/Spectral or no static solution
  else if (model->opt.discretization == ASM::Spline && pOpt.empty() && !oneD)
    pOpt[SIMoptions::GLOBAL] = "Greville point projection";

  if (model->opt.discretization < ASM::Spline && !model->opt.hdf5.empty())
  {
    IFEM::cout <<"\n ** HDF5 output is available for spline discretization only"
               <<". Deactivating...\n"<< std::endl;
    model->opt.hdf5.clear();
  }

  const char* prefix[pOpt.size()];
  if (model->opt.format >= 0 || model->opt.dumpHDF5(infile))
    for (i = 0, pit = pOpt.begin(); pit != pOpt.end(); i++, pit++)
      prefix[i] = pit->second.c_str();

  Matrix eNorm, ssol;
  Vector displ, load;
  Vectors projs(pOpt.size()), gNorm;
  std::vector<Mode> modes;
  std::vector<Mode>::const_iterator it;

  if (aSim)
    aSim->setupProjections();

  DataExporter* exporter = NULL;
  if (model->opt.dumpHDF5(infile))
  {
    IFEM::cout <<"\nWriting HDF5 file "<< model->opt.hdf5
               <<".hdf5"<< std::endl;

    // Include secondary results only if no projection has been requested.
    // The secondary results will be projected anyway, but without the
    // nodal averaging across patch boundaries in case of multiple patches.
    int results = DataExporter::PRIMARY |
                  DataExporter::NORMS   |
                  DataExporter::DISPLACEMENT;
    if (pOpt.empty() && !model->opt.pSolOnly)
      results |= DataExporter::SECONDARY;

    exporter = new DataExporter(true);
    if (staticSol)
    {
      exporter->registerField("u", "solution", DataExporter::SIM, results);
      exporter->setFieldValue("u", model, aSim ? &aSim->getSolution() : &displ);
      for (i = 0, pit = pOpt.begin(); pit != pOpt.end(); i++, pit++) {
        exporter->registerField(prefix[i], "projected", DataExporter::SIM,
                                DataExporter::SECONDARY, prefix[i]);
        exporter->setFieldValue(prefix[i], model,
                                aSim ? &aSim->getProjection(i) : &projs[i]);
      }
      exporter->setNormPrefixes(prefix);
    }
    if (model->opt.eig > 0)
    {
      exporter->registerField("eig", "eigenmode", DataExporter::SIM,
                              DataExporter::EIGENMODES);
      exporter->setFieldValue("eig", model, &modes);
    }
    exporter->registerWriter(new HDF5Writer(model->opt.hdf5,
                                            model->getProcessAdm()));
    exporter->registerWriter(new XMLWriter(model->opt.hdf5,
                                           model->getProcessAdm()));
  }

  switch (iop+model->opt.eig) {
  case 0:
  case 5:
    // Static solution: Assemble [Km] and {R}
    model->setMode(SIM::STATIC);
    model->setQuadratureRule(model->opt.nGauss[0],true,true);
    model->initSystem(model->opt.solver);
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

    if (!pOpt.empty())
      IFEM::cout << std::endl;

    if (!noError)
    {
      // Evaluate solution norms
      model->setQuadratureRule(model->opt.nGauss[1]);
      if (!model->solutionNorms(Vectors(1,displ),projs,eNorm,gNorm))
        return 4;
    }

    if (!gNorm.empty())
    {
      const Vector& norm = gNorm.front();
      if (oneD)
      {
        IFEM::cout <<"L2-norm: |u^h| = (u^h,u^h)^0.5     : "<< norm(1);
        if (norm.size() > 2)
          IFEM::cout <<"\n           |u| = (u,u)^0.5         : "<< norm(2)
                     <<"\n           |e| = (u^h-u,u^h-u)^0.5 : "<< norm(3);
      }
      else
      {
        IFEM::cout <<"Energy norm |u^h| = a(u^h,u^h)^0.5   : "<< norm(1);
        if (norm(2) != 0.0)
          IFEM::cout <<"\nExternal energy ((f,u^h)+(t,u^h)^0.5 : "<< norm(2);
        if (model->haveAnaSol() && norm.size() >= 4)
          IFEM::cout <<"\nExact norm  |u|   = a(u,u)^0.5       : "<< norm(3)
                     <<"\nExact error a(e,e)^0.5, e=u-u^h      : "<< norm(4)
                     <<"\nExact relative error (%) : "
                     << norm(4)/norm(3)*100.0;
        IFEM::cout <<"Energy norm |u^h| = a(u^h,u^h)^0.5   : "<< norm(1);
        if (norm(2) != 0.0)
          IFEM::cout <<"\nExternal energy ((f,u^h)+(t,u^h)^0.5 : "<< norm(2);
        if (model->haveAnaSol() && norm.size() >= 4)
          IFEM::cout <<"\nExact norm  |u|   = a(u,u)^0.5       : "<< norm(3)
                     <<"\nExact error a(e,e)^0.5, e=u-u^h      : "<< norm(4)
                     <<"\nExact relative error (%) : "
                     << norm(4)/norm(3)*100.0;
      }
      size_t j = 1;
      for (pit = pOpt.begin(); pit != pOpt.end() && j < gNorm.size(); pit++,j++)
      {
        IFEM::cout <<"\n\n>>> Error estimates based on "<< pit->second <<" <<<";
        IFEM::cout <<"\nEnergy norm |u^r| = a(u^r,u^r)^0.5   : "<< gNorm[j](1);
        IFEM::cout <<"\nError norm a(e,e)^0.5, e=u^r-u^h     : "<< gNorm[j](2);
        IFEM::cout <<"\n- relative error (% of |u^r|) : "
                   << gNorm[j](2)/gNorm[j](1)*100.0;
        if (j == 0) continue;

        if (model->haveAnaSol())
          IFEM::cout <<"\nExact error a(e,e)^0.5, e=u-u^r      : "<< gNorm[j](5)
                     <<"\n- relative error (% of |u|)   : "
                     << gNorm[j](5)/norm(3)*100.0
                     <<"\nEffectivity index             : "
                     << gNorm[j](2)/norm(4);

        IFEM::cout <<"\nL2-norm |s^r| =(s^r,s^r)^0.5         : "<< gNorm[j](3);
        IFEM::cout <<"\nL2-error (e,e)^0.5, e=s^r-s^h        : "<< gNorm[j](4);
        IFEM::cout <<"\n- relative error (% of |s^r|) : "
                   << gNorm[j](4)/gNorm[j](3)*100.0;
      }
      IFEM::cout << std::endl;
    }

    model->dumpResults(displ,0.0,IFEM::cout,true,6);

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
    model->setQuadratureRule(model->opt.nGauss[0],true,true);
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

    for (int iStep = 1; aSim->adaptMesh(iStep); iStep++)
      if (!aSim->solveStep(infile,iStep))
        return 5;
      else if (!aSim->writeGlv(infile,iStep,4))
        return 6;
      else if (exporter)
        exporter->dumpTimeLevel(NULL,true);

  case 100:
    break; // Model check

  default:
    // Free vibration: Assemble [Km] and [M]
    model->setMode(SIM::VIBRATION);
    model->setQuadratureRule(model->opt.nGauss[0],true,true);
    model->initSystem(model->opt.solver,2,0);
    if (!model->assembleSystem())
      return 5;

    if (!model->systemModes(modes))
      return 6;
  }

  utl::profiler->start("Postprocessing");

  if (iop != 10 && model->opt.format >= 0)
  {
    int geoBlk = 0, nBlock = 0;

    // Write VTF-file with model geometry
    if (!model->writeGlvG(geoBlk,infile))
      return 7;

    // Write boundary tractions, if any
    if (!model->writeGlvT(1,geoBlk,nBlock))
      return 8;

    // Write Dirichlet boundary conditions
    if (!model->writeGlvBC(nBlock))
      return 8;

    // Write temperature field, if specified
    const RealFunc* temp;
    const LinearElasticity* lelp;
    if ((lelp = dynamic_cast<const LinearElasticity*>(model->getProblem())))
      if ((temp = lelp->getTemperature()))
        if (!model->writeGlvF(*temp,"Temperature",1,nBlock))
          return 9;

    // Write load vector to VTF-file
    if (!model->writeGlvV(load,"Load vector",1,nBlock))
      return 10;

    // Write solution fields to VTF-file
    if (!model->writeGlvS(displ,1,nBlock))
      return 11;

    // Write projected solution fields to VTF-file
    size_t i = 0;
    int iBlk = 100;
    for (pit = pOpt.begin(); pit != pOpt.end(); pit++, i++, iBlk += 10)
      if (!model->writeGlvP(projs[i],1,nBlock,iBlk,pit->second.c_str()))
        return 12;

    // Write eigenmodes
    bool isFreq = model->opt.eig==3 || model->opt.eig==4 || model->opt.eig==6;
    for (it = modes.begin(); it != modes.end(); it++)
      if (!model->writeGlvM(*it,isFreq,nBlock))
        return 13;

    // Write element norms
    if (!model->writeGlvN(eNorm,1,nBlock,prefix))
      return 14;

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
    IFEM::cout <<"\nWriting updated g2-file "<< infile << std::endl;
    model->dumpGeometry(osg);
    if (!displ.empty())
    {
      // Write solution (control point values) to ASCII files
      std::ofstream osd(strcat(strtok(infile,"."),".dis"));
      osd.precision(18);
      IFEM::cout <<"\nWriting deformation to file "<< infile << std::endl;
      utl::LogStream log(osd);
      model->dumpPrimSol(displ,log,false);
      std::ofstream oss(strcat(strtok(infile,"."),".sol"));
      oss.precision(18);
      IFEM::cout <<"\nWriting solution to file "<< infile << std::endl;
      utl::LogStream log2(oss);
      model->dumpSolution(displ,log2);
    }
    if (!modes.empty())
    {
      // Write eigenvectors to ASCII files
      std::ofstream ose(strcat(strtok(infile,"."),".eig"));
      ose.precision(18);
      IFEM::cout <<"\nWriting eigenvectors to file "<< infile << std::endl;
      utl::LogStream log(ose);
      for (it = modes.begin(); it != modes.end(); it++)
      {
        ose <<"# Eigenvector_"<< it->eigNo <<" Eigenvalue="<< it->eigVal <<"\n";
        model->dumpPrimSol(it->eigVec,log,false);
      }
    }
  }

  utl::profiler->stop("Postprocessing");
  delete theSim;
  delete exporter;
  return 0;
}
