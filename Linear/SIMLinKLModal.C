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
#include "KirchhoffLovePlate.h"
#include "KirchhoffLoveShell.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml2.h"


SIMLinKLModal::SIMLinKLModal (std::vector<Mode>& modes, bool shell)
  : SIMLinElKL(nullptr,shell,true), SIMmodal(modes)
{
  parsed = false;
  alpha1 = alpha2 = 0.0;
}


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
  if (!this->assembleModalSystem(time,mSol,Rhs,
                                 myProblem->getIntegrationPrm(2),
                                 myProblem->getIntegrationPrm(3),
                                 alpha1,alpha2))
    return false;

  // Swap the equation systems such that the dynamic simulation driver
  // operates on the modal system
  return this->swapSystem(myEqSys,mySam);
}


const Vectors& SIMLinKLModal::expandSolution (const Vectors& mSol, bool swapBck)
{
  if (!this->expandSolution(mSol,sol))
    sol.clear();

  // Swap back to the full equation system data for postprocessing
  // and assembly of load vector for the next time step
  if (swapBck)
    this->swapSystem(myEqSys,mySam);

  return sol;
}


const Vector& SIMLinKLModal::expandedSolution (int idx) const
{
  if (idx >= 0 && idx < (int)sol.size())
    return sol[idx];

  static Vector empty;
  return empty;
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
  if (parsed)
    IFEM::cout <<"\t(skipped)"<< std::endl;
  else if (!strcasecmp(elem->Value(),"newmarksolver"))
  {
    utl::getAttribute(elem,"alpha1",alpha1);
    utl::getAttribute(elem,"alpha2",alpha2);
  }
  else
    return this->SIMLinElKL::parse(elem);

  return true;
}


bool SIMLinKLModal::preprocessB ()
{
  parsed = true;
  this->setIntegrationPrm(0,alpha1);
  this->setIntegrationPrm(1,alpha2);
  return this->SIMLinElKL::preprocessB();
}
