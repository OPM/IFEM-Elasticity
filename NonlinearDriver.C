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
#include "Elasticity.h"
#include "DataExporter.h"
#include "IFEM.h"
#include "tinyxml.h"


NonlinearDriver::NonlinearDriver (SIMbase& sim, bool linear) : NonLinSIM(sim)
{
  opt.pSolOnly = true;
  calcEn = true;
  if (linear)
    iteNorm = NONE;
}

bool NonlinearDriver::parse (char* keyWord, std::istream& is)
{
  if (!strncasecmp(keyWord,"TIME_STEPPING",13))
    return params.parse(keyWord,is);
  else if (!strncasecmp(keyWord,"NO_ENERGY",9))
    calcEn = false; // switch off energy norm calculation
  else if (!strncasecmp(keyWord,"ENERGY2",7))
    calcEn = 2; // also print the square of the global norm values
  else
    return this->NonLinSIM::parse(keyWord,is);

  return true;
}


bool NonlinearDriver::parse (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"nonlinearsolver"))
  {
    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strncasecmp(child->Value(),"noEnergy",8))
        calcEn = false; // switch off energy norm calculation
      else if (!strncasecmp(child->Value(),"energy2",7))
        calcEn = 2; // also print the square of the global norm values
      else
        params.parse(child);
  }
  else if (!strcasecmp(elem->Value(),"postprocessing"))
    if (elem->FirstChildElement("direct2nd"))
      opt.pSolOnly = false;

  return this->NonLinSIM::parse(elem);
}


bool NonlinearDriver::solutionNorms (const TimeDomain& time,
                                     double zero_tol, std::streamsize outPrec)
{
  if (msgLevel < 0 || solution.empty()) return true;

  const size_t nsd = model.getNoSpaceDim();

  size_t iMax[nsd];
  double dMax[nsd];
  double normL2 = model.solutionNorms(solution.front(),dMax,iMax);

  RealArray RF;
  bool haveReac = model.getCurrentReactions(RF,solution.front());

  Vectors gNorm;
  if (calcEn)
  {
    model.setMode(SIM::RECOVERY);
    model.setQuadratureRule(opt.nGauss[1]);
    if (!model.solutionNorms(time,solution,gNorm))
      gNorm.clear();
  }

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


/*!
  This method controls the load incrementation loop of the finite deformation
  simulation. It uses the automatic increment size adjustment of the TimeStep
  class and supports iteration cut-back in case of divergence.
*/

int NonlinearDriver::solveProblem (DataExporter* writer,
                                   utl::LogStream* oss, double dtDump,
                                   double zero_tol, std::streamsize outPrec)
{
  std::streamsize normPrec = outPrec > 3 ? outPrec : 0;

  if (dtDump <= 0.0) dtDump = params.stopTime + 1.0;
  double nextDump = params.time.t + dtDump;
  double nextSave = params.time.t + opt.dtSave;
  bool getMaxVals = opt.format >= 0 && !opt.pSolOnly;
  const Elasticity* elp = dynamic_cast<const Elasticity*>(model.getProblem());
  if (!elp) getMaxVals = false;

  int iStep = 0; // Save initial state to VTF
  if (opt.format >= 0 && params.multiSteps() && params.time.dt > 0.0)
    if (!this->saveStep(-(++iStep),params.time.t))
      return 4;

  // Initialize the linear solver
  this->initEqSystem();

  SIMoptions::ProjectionMap::const_iterator pit = opt.project.begin();
  if (pit != opt.project.end() && elp) getMaxVals = true;

  // Invoke the time-step loop
  SIM::ConvStatus stat = SIM::OK;
  while (this->advanceStep(params))
  {
    do
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
      return 5;

    if (pit != opt.project.end())
    {
      // Project the secondary results onto the spline basis
      model.setMode(SIM::RECOVERY);
      if (!model.project(proSol,solution.front(),pit->first,params.time))
        return 6;
    }

    // Print solution components at the user-defined points
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

    if (params.hasReached(nextSave))
    {
      // Save solution variables to VTF for visualization
      if (opt.format >= 0)
      {
        if (!this->saveStep(++iStep,params.time.t))
          return 7;

        // Write projected solution fields to VTF-file
        if (!model.writeGlvP(proSol,iStep,nBlock,110,pit->second.c_str(),
                             elp ? elp->getMaxVals() : nullptr))
          return 8;
      }

      // Save solution variables to HDF5
      if (writer)
        if (!writer->dumpTimeLevel(&params))
          return 9;

      nextSave = params.time.t + opt.dtSave;
      if (nextSave > params.stopTime)
        nextSave = params.stopTime; // Always save the final step
    }
    else if (getMaxVals)
    {
      if (!model.eval2ndSolution(solution.front(),params.time.t))
        return 10;

      if (!model.evalProjSolution(proSol,*elp->getMaxVals()))
        return 11;
    }

    // Print out the maximum von Mises stress, etc., if present
    if (getMaxVals && myPid == 0)
    {
      size_t id = model.getNoSpaceDim()*2 + 1;
      elp->printMaxVals(outPrec,id);   // von Mises stress
      elp->printMaxVals(outPrec,id+1); // plastic strain, Epp
      elp->printMaxVals(outPrec,id+6); // stress triaxiality, T
      elp->printMaxVals(outPrec,id+7); // Lode parameter, L
    }
  }

  return 0;
}
