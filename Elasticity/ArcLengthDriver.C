// $Id$
//==============================================================================
//!
//! \file ArcLengthDriver.C
//!
//! \date Apr 12 2018
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Arc-length solution driver for nonlinear isogeometric FEM simulators.
//!
//==============================================================================

#include "ArcLengthDriver.h"
#include "SIMoutput.h"
#include "IFEM.h"
#include "IntegrandBase.h"
#include "Profiler.h"
#include "Utilities.h"
#include "tinyxml2.h"


ArcLengthDriver::ArcLengthDriver (SIMbase& sim, bool, bool adaptive)
  : NonlinearDriver(sim,false,adaptive)
{
  nRHSvec = 2; // 0: residual force, 1: load gradient

  // Default arc-length parameters
  arclen = 0.1;
  beta   = -1.0;
}


bool ArcLengthDriver::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),"nonlinearsolver"))
  {
    const char* value;
    const tinyxml2::XMLElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if ((value = utl::getValue(child,"beta")))
        beta = atof(value);
      else if ((value = utl::getValue(child,"arclen")))
        arclen = atof(value);
  }

  return this->NonlinearDriver::parse(elem);
}


bool ArcLengthDriver::advanceStep (TimeStep& param, bool)
{
  if (param.finished())
  {
    IFEM::cout <<"\n  Path integration completed."<< std::endl;
    return false;
  }

  ++param.step;
  return this->NonlinearDriver::advanceStep(param,false);
}


SIM::ConvStatus ArcLengthDriver::solveStep (TimeStep& param,
                                            SIM::SolutionMode,
                                            double zero_tolerance,
                                            std::streamsize outPrec)
{
  PROFILE1("ArcLengthDriver::solveStep");

#ifdef INT_DEBUG
  std::cout <<"\n=== ArcLengthDriver::solveStep("<< param.step
            <<"): lambda = "<< param.time.t << std::endl;
#endif

  if (solution.empty() || !model.setMode(SIM::ARCLEN))
    return SIM::FAILURE;

  // Assemble new tangent stiffness matrix, K^(0)
  // the residual forces, R^(0) and load gradient, dF/dlam
  param.iter = 0; // for the predictor step (iteration 0)
  model.setQuadratureRule(opt.nGauss[0],true);
  if (!model.assembleSystem(param.time,solution))
    return model.getProblem()->diverged() ? SIM::DIVERGED : SIM::FAILURE;

  // Solve for residual and tangent displacements
  double lgNorm = this->solveLinearizedSystem(param.time.t);
  if (lgNorm < 0.0)
    return SIM::FAILURE;

  // Calculate the incremental Load Proportionality Factor of the predictor step
  double dellam, tanNorm = tansol.norm2();
  if (beta < 0.0)
    dellam = arclen / sqrt(1.0 + tanNorm*tanNorm);
  else if (beta > 0.0)
    dellam = arclen / sqrt(beta*lgNorm*lgNorm + tanNorm*tanNorm);
  else
    dellam = arclen / tanNorm;

  if (!incsol.empty() && incsol.dot(tansol) < 0.0)
    dellam = -dellam;

  // Update the total Load Proportionality Factor
  param.time.t += dellam;
  param.time.dt = dellam;

  // Calculate the incremental solution vector of the predictor step
  linsol.add(tansol,dellam);
  incsol = linsol;
  // Update the total solution
  if (!this->updateConfiguration(param))
    return SIM::FAILURE;

  // Assemble new tangent stiffness matrix, K^(0)
  // the residual forces, R^(0) and load gradient, dF/dlam
  param.time.first = false; // for the updated configuration
  if (!model.assembleSystem(param.time,solution))
    return model.getProblem()->diverged() ? SIM::DIVERGED : SIM::FAILURE;

  if (msgLevel >= 0)
    model.printStep(param.step,param.time);

  // Extract the residual force vector of the predictor step
  if (!model.extractLoadVec(residual,0))
    return SIM::FAILURE;

  // Newton-Raphson iteration loop
  bool poorConvg = false;
  while (param.iter <= maxit)
    switch (this->checkConvergence(param))
      {
      case SIM::CONVERGED:
        if (!this->solutionNorms(param.time,zero_tolerance,outPrec))
          return SIM::FAILURE;
        return SIM::CONVERGED;

      case SIM::DIVERGED:
        return SIM::DIVERGED;

      case SIM::SLOW:
        poorConvg = true;

      default:
        // Not convereged yet, solve for new residual and tangent displacements
        if (this->solveLinearizedSystem(param.time.t,++param.iter) < 0.0)
          return SIM::FAILURE;

        // Calculate the iterative Load Proportionality Factor
        if (beta < 0.0)
          dellam = -incsol.dot(linsol) / (incsol.dot(tansol) + param.time.dt);
        else
          dellam = -tansol.dot(linsol) / (1.0 + tansol.dot(tansol));

        // Update the total and incremental Load Proportionality Factor
        param.time.t  += dellam;
        param.time.dt += dellam;

        // Update the incremental solution vector
        linsol.add(tansol,dellam);
        incsol.add(linsol);
        // Update the total solution
        if (!this->updateConfiguration(param))
          return SIM::FAILURE;

        // Assemble new tangent stiffness matrix, K^(i) (if full NR)
        // the residual forces, R^(i) and updated load gradient, dF/dlam
        bool newTangent = param.iter <= nupdat;
        model.setMode(param.iter > nupdat ? SIM::RHS_ONLY : SIM::ARCLEN);
        if (!model.assembleSystem(param.time,solution,newTangent,poorConvg))
          return model.getProblem()->diverged() ? SIM::DIVERGED : SIM::FAILURE;

        poorConvg = false;
      }

  return SIM::DIVERGED;
}


double ArcLengthDriver::solveLinearizedSystem (double lambda, int iter)
{
#ifdef INT_DEBUG
  std::cout <<"\n=== ArcLengthDriver::solveLinearizedSystem("<< iter
            <<"): lambda = "<< lambda << std::endl;
#endif

  if (!model.extractLoadVec(residual,0))
    return -0.1;

  double fgNorm = 0.0;
  if (beta > 0.0)
  {
    Vector FextGrad;
    if (model.extractLoadVec(FextGrad,1))
      fgNorm = FextGrad.norm2();
  }

  double* rCondPtr = rCond < 0.0 ? nullptr : &rCond;
  if (!model.solveSystem(linsol,msgLevel-1,rCondPtr,"residual disp",0))
    return -1.0;

  if (!model.solveSystem(tansol,msgLevel-1,nullptr,"tangential disp",1))
    return -2.0;

  return fgNorm;
}


SIM::ConvStatus ArcLengthDriver::checkConvergence (TimeStep& param)
{
  static double convTol   = 0.0;
  static double prevNorm  = 0.0;
  static int    nIncrease = 0;

  SIM::ConvStatus status = SIM::OK;
  double enorm, resNorm, linsolNorm;
  model.iterationNorms(linsol,residual,enorm,resNorm,linsolNorm);
  double norm = iteNorm == ENERGY ? enorm : resNorm;
  if (iteNorm == L2SOL) norm = linsolNorm;

  if (param.iter == 0)
  {
    if (linsolNorm == 0.0)
      return SIM::CONVERGED; // No load on this step

    if (refNopt == ALL || fabs(norm) > refNorm)
      refNorm = fabs(norm);

    if (refNorm*rTol > aTol) {
      convTol = rTol;
      prevNorm = (norm /= refNorm);
    }
    else {
      convTol = aTol;
      refNorm = 1.0;
      prevNorm = norm;
    }

    nIncrease = 0;
  }
  else
    norm /= refNorm;

  // Check for slow convergence
  if (param.iter > 1 && prevNorm > 0.0 && fabs(norm) > prevNorm*0.1)
    status = SIM::SLOW;

  if (msgLevel > 0)
  {
    // Print convergence history
    utl::LogStream& cout = model.getProcessAdm().cout;
    std::ios::fmtflags oldFlags = cout.flags(std::ios::scientific);
    std::streamsize oldPrec = cout.precision(3);
    cout <<"  iter="<< param.iter
         <<"  time="<< param.time.t
         <<"  conv="<< fabs(norm)
         <<"  enen="<< enorm
         <<"  resn="<< resNorm
         <<"  incn="<< linsolNorm;
    if (rCond > 0.0)
      cout <<"  cond="<< 1.0/rCond;
    cout << std::endl;

    // Find and print out the worst DOF(s) when detecting slow convergence
    if (status == SIM::SLOW && prnSlow > 0)
      this->printWorst(cout,convTol*refNorm);

    cout.flags(oldFlags);
    cout.precision(oldPrec);
  }

  // Check for convergence or divergence
  if (fabs(norm) < convTol && (param.iter > 0 || refNopt == ALL))
    status = SIM::CONVERGED;
  else if (std::isnan(linsolNorm))
    status = SIM::DIVERGED;
  else if (fabs(norm) <= fabs(prevNorm))
    nIncrease = 0;
  else if (++nIncrease > maxIncr || fabs(norm) > divgLim)
    status = SIM::DIVERGED;

  prevNorm = norm;
  return status;
}
