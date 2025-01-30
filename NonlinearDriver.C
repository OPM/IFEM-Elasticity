// $Id$
//==============================================================================
//!
//! \file NonlinearDriver.C
//!
//! \date Jul 15 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Nonlinear driver for isogeometric finite deformation FEM analysis.
//!
//==============================================================================

#include "NonlinearDriver.h"
#include "SIMoutput.h"
#include "AdaptiveSetup.h"
#include "ASMunstruct.h"
#include "Elasticity.h"
#include "DataExporter.h"
#include "HDF5Restart.h"
#include "Profiler.h"
#include "IFEM.h"
#include "tinyxml2.h"


NonlinearDriver::NonlinearDriver (SIMbase& sim, bool linear, bool adaptive)
  : NonLinSIM(sim, linear ? NONE : ENERGY)
{
  aStep = 0;
  save0 = opt.pSolOnly = true;

  if (adaptive)
  {
    adap = new AdaptiveSetup(static_cast<SIMoutput&>(sim));
    calcEn = 0;
  }
  else
  {
    adap = nullptr;
    calcEn = 1;
  }
}


NonlinearDriver::~NonlinearDriver ()
{
  delete adap;
}


bool NonlinearDriver::read (const char* fileName)
{
  if (adap) inpfile = fileName;
  return this->NonLinSIM::read(fileName);
}


bool NonlinearDriver::parse (char* keyWord, std::istream& is)
{
  if (!strncasecmp(keyWord,"TIME_STEPPING",13))
    return params.parse(keyWord,is);
  else if (!strncasecmp(keyWord,"NO_ENERGY",9))
    calcEn = 0; // switch off energy norm calculation
  else if (!strncasecmp(keyWord,"ENERGY2",7))
    calcEn = 2; // also print the square of the global norm values
  else if (!strncasecmp(keyWord,"ADAPTIVE",8) && adap)
    return adap->parse(keyWord,is);
  else
    return this->NonLinSIM::parse(keyWord,is);

  return true;
}


bool NonlinearDriver::parse (const tinyxml2::XMLElement* elem)
{
  if (adap && !strcasecmp(elem->Value(),"adaptive"))
    return adap->parse(elem);

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  if (!strcasecmp(elem->Value(),"nonlinearsolver"))
  {
    for (; child; child = child->NextSiblingElement())
      if (!strncasecmp(child->Value(),"noEnergy",8))
        calcEn = 0; // switch off energy norm calculation
      else if (!strncasecmp(child->Value(),"energy2",7))
        calcEn = 2; // also print the square of the global norm values
      else
        params.parse(child);
  }

  else if (!strcasecmp(elem->Value(),"postprocessing"))
  {
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"direct2nd"))
        opt.pSolOnly = false;
      else if (!strcasecmp(child->Value(),"skipInit"))
        save0 = false;
      else if (!strcasecmp(child->Value(),"resultpoints") &&
               child->FirstChildElement("grid"))
        save0 = false; // Deactivate initial configuration dump if grid output
  }

  return params.parse(elem) && this->NonLinSIM::parse(elem);
}


bool NonlinearDriver::solutionNorms (const TimeDomain& time,
                                     double zero_tol, std::streamsize outPrec)
{
  if (msgLevel < 0 || solution.empty()) return true;

  const size_t nsd = model.getNoSpaceDim();

  size_t iMax[3];
  double dMax[3];
  RealArray RF, Fext;
  double normL2 = model.solutionNorms(solution.front(),dMax,iMax);
  bool haveReac = model.getCurrentReactions(RF,solution.front());

  if (calcEn)
  {
    model.setMode(SIM::NORMS);
    model.setQuadratureRule(opt.nGauss[1]);
    if (!model.solutionNorms(time,solution,gNorm))
      gNorm.clear();
  }
  else
    gNorm.clear();

  if (myPid > 0) return true;

  std::streamsize stdPrec = outPrec > 0 ? IFEM::cout.precision(outPrec) : 0;
  double old_tol = utl::zero_print_tol;
  utl::zero_print_tol = zero_tol;

  IFEM::cout <<"  Primary solution summary: L2-norm            : "
             << utl::trunc(normL2);

  for (unsigned char d = 0; d < nsd; d++)
    if (utl::trunc(dMax[d]) != 0.0)
      IFEM::cout <<"\n                            Max "<< char('X'+d)
                 <<"-displacement : "<< dMax[d] <<" node "<< iMax[d];

  if (model.extractScalars(Fext))
  {
    IFEM::cout <<"\n  Total external load: Sum(Fex) =";
    for (double f : Fext) IFEM::cout <<" "<< utl::trunc(f);
  }

  if (haveReac)
  {
    IFEM::cout <<"\n  Total reaction forces: Sum(R) =";
    for (size_t i = 1; i < RF.size(); i++)
      IFEM::cout <<" "<< utl::trunc(RF[i]);
    if (utl::trunc(RF.front()) != 0.0)
      IFEM::cout <<"\n  displacement*reactions: (R,u) = "<< RF.front();
  }

  if (!gNorm.empty())
    this->printNorms(gNorm.front(),IFEM::cout);

  IFEM::cout << std::endl;
  utl::zero_print_tol = old_tol;
  if (stdPrec > 0) IFEM::cout.precision(stdPrec);
  return true;
}


void NonlinearDriver::printNorms (const Vector& norm, utl::LogStream& os) const
{
  if (norm.size() > 0)
  {
    os <<"\n  Energy norm:    |u^h| = a(u^h,u^h)^0.5 : "<< utl::trunc(norm(1));
    if (calcEn == 2)
    {
      std::streamsize oldPrec = os.precision(10);
      os <<"\t a(u^h,u^h) = "<< utl::trunc(norm(1)*norm(1));
      os.precision(oldPrec);
    }
  }
  if (norm.size() > 1 && utl::trunc(norm(2)) != 0.0)
  {
    os <<"\n  External energy: ((f,u^h)+(t,u^h))^0.5 : "<< norm(2);
    if (calcEn == 2)
    {
      std::streamsize oldPrec = os.precision(10);
      os <<"\t(f,u)+(t,u) = "<< norm(2)*norm(2);
      os.precision(oldPrec);
    }
  }
  if (norm.size() > 2)
    os <<"\n  Stress norm, L2: (sigma^h,sigma^h)^0.5 : "<< norm(3);
  if (norm.size() > 3)
    os <<"\n  Pressure norm, L2:       (p^h,p^h)^0.5 : "<< norm(4)
       <<"\t(p^h = trace(sigma^h)/3)";
  if (norm.size() > 4)
    os <<"\n  Deviatoric stress norm:  (s^d,s^d)^0.5 : "<< norm(5)
       <<"\t(s^d = sigma^h - p^h*I)";
  if (norm.size() > 5)
    os <<"\n  Stress norm, von Mises: vm(sigma^h)    : "<< norm(6);
}


bool NonlinearDriver::calcInterfaceForces (double t)
{
  bool ok = model.assembleForces(solution.front(),t,&myReacts,&myForces);
  if (ok)
    model.printIFforces(myForces,myWeights);
  else
    std::cerr <<" *** Interface force calculation failes."<< std::endl;

  return ok;
}


/*!
  This method controls the load incrementation loop of the nonlinear simulation.
  It uses the automatic increment size adjustment of the TimeStep class
  and supports iteration cut-back in case of divergence.
  It also supports adaptive mesh refinement based on error indicators.
*/

int NonlinearDriver::solveProblem (DataExporter* writer, HDF5Restart* restart,
                                   utl::LogStream* oss, bool printMax,
                                   double dtDump, double zero_tol,
                                   std::streamsize outPrec)
{
  std::streamsize normPrec = outPrec > 3 ? outPrec : 0;

  SIMoptions::ProjectionMap::const_iterator pit = opt.project.begin();
  bool doProject = pit != opt.project.end();
  if (doProject) proSol.resize(1);

  if (dtDump <= 0.0) dtDump = params.stopTime + 1.0;
  double nextDump = params.time.t + dtDump;
  double nextSave = params.time.t + opt.dtSave;
  bool getMaxVals = opt.format >= 0 && !opt.pSolOnly;
  const Elasticity* elp = dynamic_cast<const Elasticity*>(model.getProblem());
  if (!elp)
    getMaxVals = printMax = false;
  else if (doProject)
    getMaxVals = true;

  int iStep = aStep = 0; // Save initial state to VTF
  if (save0 && opt.format >= 0 && params.multiSteps() && params.time.dt > 0.0)
    if (!this->saveStep(-(++iStep),params.time.t))
      return 4;

  // Initialize mesh adaptation parameters
  if (adap && !adap->initPrm(1))
    return 5;

  // Start the load incrementation loop
  SIM::ConvStatus stat = SIM::OK;
  while (this->advanceStep(params))
  {
    int bStep = aStep; // Check for mesh adaptation
    if (params.step > 1 && !this->adaptMesh(aStep))
      return 6;
    else if (aStep > bStep)
      IFEM::cout <<"\nResuming nonlinear solution on the new mesh"<< std::endl;

    do // Cut-back loop
    {
      if (stat == SIM::DIVERGED)
      {
        // Try cut-back with a smaller time step when diverging
        if (!params.cutback()) break;

        std::copy(solution[1].begin(),solution[1].end(),solution[0].begin());
        model.updateConfiguration(solution.front());
        refNorm = 1.0; // Reset the reference norm
      }

      // Solve the nonlinear FE problem at this load step
      stat = this->solveStep(params,SIM::STATIC,zero_tol,normPrec);
    }
    while (stat == SIM::DIVERGED);

    if (stat != SIM::CONVERGED)
      return 7;

    if (model.haveBoundaryReactions())
      if (!this->calcInterfaceForces(params.time.t))
        return 9;

    if (doProject)
    {
      // Project the secondary results onto the spline basis
      model.setMode(SIM::RECOVERY);
      Matrix projs(proSol.front());
      if (!model.project(projs,solution.front(),pit->first,params.time))
        return 8;
    }

    if (adap)
    {
      // Evaluate error norms
      model.setMode(SIM::NORMS);
      model.setQuadratureRule(opt.nGauss[1]);
      if (!model.solutionNorms(params.time,solution,proSol,gNorm,&eNorm))
        return 9;

      IFEM::cout << std::endl;
      adap->printNorms(gNorm,Vectors(),eNorm);
    }

    utl::profiler->start("Postprocessing");

    // Print solution components at the user-defined points
    model.setMode(SIM::RECOVERY);
    this->dumpResults(params.time.t,IFEM::cout,outPrec);

    if (params.hasReached(nextDump))
    {
      // Dump primary solution for inspection or external processing
      if (oss)
        this->dumpStep(params.step,params.time.t,*oss,false);
      else
        this->dumpStep(params.step,params.time.t,IFEM::cout);

      nextDump = params.time.t + dtDump;
    }

    if (opt.dtSave <= 0.0 || params.hasReached(nextSave))
    {
      ++iStep;

      // Save solution variables to VTF for visualization
      if (opt.format >= 0)
      {
        if (!this->saveStep(iStep,params.time.t))
          return 11;

        if (!myForces.empty())
          if (!model.writeGlvV(myForces,"Internal forces",iStep,nBlock,2))
            return 11;

        if (doProject)
        {
          // Write projected solution fields to VTF-file
          if (!model.writeGlvP(proSol.front(),iStep,
                               nBlock,110,pit->second.c_str(),
                               elp ? elp->getMaxVals() : nullptr))
            return 11;

          // Write element norms
          if (!model.writeGlvN(eNorm,iStep,nBlock,{pit->second}))
            return 11;
        }
      }

      if (elp) elp->enableMaxValCalc(false);

      // Save solution variables to HDF5 file
      if (writer && !writer->dumpTimeLevel(&params))
        return 12;

      // Save solution state to restart HDF5 file
      if (restart && restart->dumpStep(params))
      {
        SerializeMap data;
        if (this->serialize(data) && !restart->writeData(data))
          return 12;
      }

      // Save solution variables to grid files, if specified
      if (!model.saveResults(solution,params.time.t,iStep))
        return 13;

      if (elp) elp->enableMaxValCalc(true);

      nextSave = params.time.t + opt.dtSave;
      if (nextSave > params.stopTime)
        nextSave = params.stopTime; // Always save the final step
    }
    else if (getMaxVals)
    {
      if (!model.eval2ndSolution(solution.front(),params.time.t))
        return 14;

      if (!model.writeGlvP(proSol.front(),0,nBlock,0,nullptr,elp->getMaxVals()))
        return 15;
    }

    // Print out the maximum von Mises stress, etc., if present
    if (printMax && myPid == 0)
      elp->printMaxVals(outPrec); // Print all components
    else if (getMaxVals && myPid == 0)
    {
      size_t id = elp->getNoFields(3) + 1;
      elp->printMaxVals(outPrec,id);   // von Mises stress
      elp->printMaxVals(outPrec,id+1); // plastic strain, Epp
      elp->printMaxVals(outPrec,id+6); // stress triaxiality, T
      elp->printMaxVals(outPrec,id+7); // Lode parameter, L
    }

    utl::profiler->stop("Postprocessing");
  }

  return 0;
}


bool NonlinearDriver::serialize (SerializeMap& data) const
{
  return params.serialize(data) && this->NonLinSIM::serialize(data);
}


bool NonlinearDriver::deSerialize (const SerializeMap& data)
{
  return params.deSerialize(data) && this->NonLinSIM::deSerialize(data);
}


bool NonlinearDriver::adaptMesh (int& aStep)
{
  if (!adap) return true; // No mesh-refinement, silently ignore

  // Check for adaptive mesh refinement and extract refinement indicators
  LR::RefineData prm;
  int ierr = adap->calcRefinement(prm,aStep+1,gNorm,eNorm.getRow(adap->eIdx()));
  if (ierr < 0)
    return false;
  else if (ierr == 0)
    return true; // The mesh is fine, continue simulation on current mesh

  // Save the size of the solution vector array for the solution transfer log,
  // because refine() will resize it according to the refined mesh
  size_t nsol = solution.size();
  size_t nsv1 = solution.empty() ? 0 : solution.front().size();

  // Do the mesh refinement
  bool ok = model.refine(prm,solution);
  // Write mesh files for inspection, if requested
  adap->writeMesh(++aStep);

  // Read the input file again to set up the refined model
  model.clearProperties();
  if (!ok || !model.readModel(inpfile.c_str()))
    return false;

  if (!model.preprocess())
    return false;

  if (!this->initEqSystem(true,model.getNoFields()))
    return false;

  // Transfer primary solution variables onto the new mesh
  IFEM::cout <<"\nTransferring ";
  if (nsol > 1) IFEM::cout << nsol <<"x";
  IFEM::cout << nsv1 <<" solution variables to the new mesh"<< std::endl;
  Vectors soli(nsol,Vector(model.getNoDOFs()));
  for (size_t i = 0; i < nsol; i++)
    for (int p = 0; p < model.getNoPatches(); p++)
      model.injectPatchSolution(soli[i],solution[p*nsol+i],model.getPatch(p+1));

  // Write updated geometry and (homogeneous) Dirichlet BCs to VTF-file
  return opt.format < 0 ? true : this->saveModel(geoBlk,nBlock);
}
