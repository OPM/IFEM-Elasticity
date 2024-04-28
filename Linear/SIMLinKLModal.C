// $Id$
//==============================================================================
//!
//! \file SIMLinKLModal.C
//!
//! \date Nov 21 2023
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for modal linear Kirchhoff-Love FEM analysis.
//!
//==============================================================================

#include "SIMLinKLModal.h"
#include "IntegrandBase.h"


bool SIMLinKLModal::assembleSystem (const TimeDomain& time,
                                    const Vectors& mSol, bool, bool)
{
  // Assemble the eigenvalue system
  if (myProblem->getMode() == SIM::VIBRATION)
    return this->SIM2D::assembleSystem(time,Vectors());

  if (time.it > 0)
    // Swap back to the full equation system for assembly of load vector
    this->swapSystem(myEqSys,mySam);

  else
  {
    // Assemble the load vector of this time step.
    // We need to do this in the first iteration only, as for linear systems
    // the load vector is not supposed to change during the iterations.
    if (!this->SIM2D::assembleSystem(time,sol,false))
      return false;

    // Extract the load vector in DOF-order
    if (!this->extractLoadVec(Rhs))
      return false;
  }

  // Assemble the modal equation system
  if (!this->assembleModalSystem(time,mSol,
                                 myProblem->getIntegrationPrm(2),
                                 myProblem->getIntegrationPrm(3)))
    return false;

  // Swap the equation systems such that the dynamic simulation driver
  // operates on the modal system
  return this->swapSystem(myEqSys,mySam);
}


const Vectors& SIMLinKLModal::expandSolution (const Vectors& mSol, bool swapBck)
{
  // Swap back to the full equation system data for postprocessing
  // and assembly of load vector for the next time step
  if (swapBck)
    this->swapSystem(myEqSys,mySam);

  return this->expandSolution(mSol);
}


bool SIMLinKLModal::serialize (std::map<std::string,std::string>& data) const
{
  return this->saveModes(data);
}


bool SIMLinKLModal::deSerialize (const std::map<std::string,std::string>& data)
{
  return this->restoreModes(data);
}


bool SIMLinKLModal::projectModes (Matrices& sesol,
                                  std::vector<std::string>& names,
                                  SIMoptions::ProjectionMethod pMethod)
{
  sesol.resize(myModes.size());
  names.resize(myProblem->getNoFields(2));
  for (size_t c = 0; c < names.size(); c++)
    names[c] = myProblem->getField2Name(c);

  bool ok = this->setMode(SIM::RECOVERY);
  for (size_t i = 0; i < myModes.size() && ok; i++)
    ok = this->project(sesol[i],myModes[i].eigVec,pMethod);

  return ok;
}


bool SIMLinKLModal::parse (const tinyxml2::XMLElement* elem)
{
  return this->parseParams(elem) || this->SIMLinElKL::parse(elem);
}


bool SIMLinKLModal::preprocessB ()
{
  parsed = true;
  this->setIntegrationPrm(0,alpha1);
  this->setIntegrationPrm(1,alpha2);
  return this->SIMLinElKL::preprocessB();
}
