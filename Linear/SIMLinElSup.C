// $Id$
//==============================================================================
//!
//! \file SIMLinElSup.C
//!
//! \date Jun 12 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for linear elastic superelement FEM analysis.
//!
//==============================================================================

#include "SIMLinElSup.h"
#include "SIMLinEl.h"
#include "SAM.h"
#include "ElementBlock.h"
#include "Utilities.h"
#include "Tensor.h"
#include "VTF.h"


SIMLinElSup::~SIMLinElSup()
{
  for (std::pair<const std::string,FEmodel>& sub : mySubSim)
  {
    delete sub.second.sim;
    delete sub.second.blk;
  }
}


bool SIMLinElSup::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"elasticity"))
    return this->SIMsupel::parse(elem);

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(),"superelement"))
    {
      std::string supId;
      if (utl::getAttribute(child,"id",supId) && !supId.empty())
      {
        // This <elasticity> block defines the FE model of a superelement
        SIMgeneric* sup = new SIMLinEl3D("Linear elastic superelement",false);
        if (sup->parse(elem))
          mySubSim[supId].sim = sup;
        else
          delete sup;
      }
      break;
    }

  return true;
}


bool SIMLinElSup::preprocessB ()
{
  // Set up links to the underlying FE model for each superelement
  for (SuperElm& sup : mySups)
  {
    std::map<std::string,FEmodel>::const_iterator sit = mySubSim.find(sup.id);
    if (sit != mySubSim.end()) sup.sim = sit->second.sim;
  }

  // Preprocess the superelement FE models for the recovery process
  bool ok = true;
  for (std::pair<const std::string,FEmodel>& sub : mySubSim)
    ok &= (sub.second.sim->preprocess({},fixDup) &&
           sub.second.sim->initSystem(LinAlg::SPARSE,0,0) &&
           sub.second.sim->setMode(SIM::RECOVERY));

  return ok;
}


ElementBlock* SIMLinElSup::tesselatePatch (size_t pidx) const
{
  if (pidx >= mySups.size() || !mySups[pidx].sim)
    return this->SIMsupel::tesselatePatch(pidx); // use simplified tesselation

  ElementBlock* supblk = nullptr;
  FEmodel& sub = const_cast<SIMLinElSup*>(this)->mySubSim[mySups[pidx].id];

  if (!sub.blk && sub.sim)
    // Create an ElementBlock consisting of all patches in the superelement
    for (ASMbase* pch : sub.sim->getFEModel())
      if (pch && !pch->empty())
      {
        ElementBlock newblk(8);
        if (!pch->tesselate(newblk,sub.sim->opt.nViz))
        {
          delete sub.blk;
          sub.blk = nullptr;
          break;
        }
        else if (!sub.blk)
          sub.blk = new ElementBlock(newblk);
        else // Append to sub.blk, without checking for unique nodal point
          sub.blk->merge(newblk,false);
      }

  if (sub.blk)
  {
    // Make a copy of sub.blk and apply the superelement transformation to it
    supblk = new ElementBlock(*sub.blk);
    supblk->transform(mySups[pidx].MVP);
  }

  return supblk;
}


bool SIMLinElSup::recoverInternalDispl (const Vector& glbSol)
{
  size_t pidx = 0;
  for (SuperElm& sup : mySups)
  {
    ASMbase* pch = myModel[pidx++];

    // Extract superelement solution vector from the global solution vector
    Vector supSol;
    pch->extractNodalVec(glbSol, supSol, mySam->getMADOF());
#if INT_DEBUG > 2
    std::cout <<"\nSolution vector for superelement "<< pch->idx+1 << supSol;
#endif

    if (sup.sim)
    {
      bool ok = true;
      Vector sol(supSol);
      if (!sup.MVP.empty()) // Transform to local superelement axes
        for (size_t i = 1; i < sol.size() && ok; i += 3)
          ok = utl::transform(sol,sup.MVP,i,true);

      // Recover the internal displacement state
      ok &= sup.sim->recoverInternals(sol,sup.sol);

      if (!sup.MVP.empty()) // Transform back to global axes
        for (size_t i = 1; i < sup.sol.size() && ok; i += 3)
          ok = utl::transform(sup.sol,sup.MVP,i);

      if (!ok)
      {
        std::cerr <<"\n *** SIMLinElSup::recoverInternalDispl: Failed to"
                  <<" recover internal displacements for superelement "
                  << pch->idx+1 << std::endl;
        return false;
      }
    }
    else // No substructure FE model - just use the superelement displacements
      sup.sol = supSol;
  }

  return true;
}


IntegrandBase* SIMLinElSup::getMyProblem () const
{
  for (const std::pair<const std::string,FEmodel>& sub : mySubSim)
    if (sub.second.sim)
      return const_cast<IntegrandBase*>(sub.second.sim->getProblem());

  return myProblem;
}


/*!
  \brief Static helper to write out scalar fields to VTF-file.
*/

static bool writeFields (const Matrix& field, int geomID,
                         int& nBlock, std::vector<IntVec>& sID, VTF* vtf)
{
  for (size_t j = 1; j <= field.rows(); j++)
    if (!vtf->writeNres(field.getRow(j), ++nBlock, geomID))
      return false;
    else if (j <= sID.size())
      sID[j-1].push_back(nBlock);
    else
      sID.push_back({nBlock});

  return true;
}


int SIMLinElSup::writeGlvS1 (const Vector& psol, int iStep, int& nBlock,
                             double, const char*, int idBlock, int, bool)
{
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return 0;
  else if (psol.empty())
    return 0;

  // Recover internal displacements for all superelements

  if (!this->recoverInternalDispl(psol))
    return -9;

  VTF* vtf = this->getVTF();
  if (!vtf) return -99;

  IntVec vID;
  std::vector<IntVec> sID;
  sID.reserve(this->getNoFields());

  int geomID = this->getStartGeo();
  for (size_t pidx = 0; pidx < myModel.size(); pidx++)
  {
    Matrix field, subfield;

    if (mySups[pidx].sim)
      for (ASMbase* pch : mySups[pidx].sim->getFEModel())
      {
        // Extract displacement vector for this sub-patch
        Vector pchvec;
        pch->extractNodalVec(mySups[pidx].sol, pchvec,
                             mySups[pidx].sim->getSAM()->getMADOF());
#if INT_DEBUG > 2
        std::cout <<"\nSolution vector for sub-patch "<< pch->idx+1 << pchvec;
#endif

        // Evaluate the internal displacement field on this sub-patch
        if (!pch->evalSolution(subfield, pchvec, mySups[pidx].sim->opt.nViz))
          return -5;

        if (field.empty())
          field = subfield;
        else
          field.augmentCols(subfield);
      }

    else // Evaluate displacement field on supernodes only
      if (!myModel[pidx]->evalSolution(field, mySups[pidx].sol, opt.nViz))
        return -1;

    if (msgLevel > 1)
      IFEM::cout <<"Writing primary solution for patch "
                 << myModel[pidx]->idx+1 <<" ("<< field.rows()
                 <<","<< field.cols() <<")"<< std::endl;

    // Output as vector field
    if (!vtf->writeVres(field, ++nBlock, ++geomID, this->getNoSpaceDim()))
      return -2;
    else
      vID.push_back(nBlock);

    // Output as scalar fields
    if (!writeFields(field, geomID, nBlock, sID, vtf))
      return -3;
  }

  // Write result block identifications

  bool ok = vID.empty() || vtf->writeDblk(vID,"Displacement",idBlock,iStep);
  for (size_t i = 0; i < sID.size() && !sID[i].empty() && ok; i++)
    ok = vtf->writeSblk(sID[i], this->getMyProblem()->getField1Name(i).c_str(),
                        idBlock++, iStep);

  return ok ? idBlock : -4;
}


/*!
  This method assumes that the internal displacements on the superelements
  already have been recovered, by calling writeGlvS1() first. Therefore,
  it does not need access to the global primary solution vector.
*/

bool SIMLinElSup::writeGlvS2 (const Vector&, int iStep, int& nBlock,
                              double, int idBlock, int)
{
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return true;

  VTF* vtf = this->getVTF();
  if (!vtf) return false;

  IntegrandBase* problem = this->getMyProblem();

  std::vector<IntVec> sID;
  sID.reserve(problem->getNoFields(2));

  int geomID = this->getStartGeo();
  for (size_t pidx = 0; pidx < myModel.size(); pidx++)
  {
    SIMoutput* supsim = mySups[pidx].sim;
    Matrix field, subfield;

    if (supsim)
      for (ASMbase* pch : supsim->getFEModel())
      {
        // Extract displacement vector for this sub-patch
        pch->extractNodalVec(mySups[pidx].sol, problem->getSolution(),
                             supsim->getSAM()->getMADOF());

        // Direct evaluation of secondary solution variables
        LocalSystem::patch = pch->idx+1;
        supsim->setPatchMaterial(pch->idx);
        if (!pch->evalSolution(subfield, *problem, supsim->opt.nViz))
          return false;

        if (field.empty())
          field = subfield;
        else
          field.augmentCols(subfield);
      }
    else // Pad with zero values for the supernodes
      field.resize(problem->getNoFields(2),
                   myModel[pidx]->getNoNodes()+1);

    if (msgLevel > 1)
      IFEM::cout <<"Writing secondary solution for patch "
                 << myModel[pidx]->idx+1 <<" ("<< field.rows()
                 <<","<< field.cols() <<")"<< std::endl;

    // Output as scalar fields
    if (!writeFields(field, ++geomID, nBlock, sID, vtf))
      return false;
  }

  // Write result block identifications

  if (idBlock <= 20)
    idBlock = 21; // avoid conflict with blocks from the 1D simulator
  bool ok = true;
  for (size_t i = 0; i < sID.size() && !sID[i].empty() && ok; i++)
    ok = vtf->writeSblk(sID[i], problem->getField2Name(i).c_str(),
                        idBlock++, iStep);

  return ok;
}
