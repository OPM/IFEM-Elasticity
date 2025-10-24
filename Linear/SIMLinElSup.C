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


bool SIMLinElSup::parse (const tinyxml2::XMLElement* elem)
{
  if (strcasecmp(elem->Value(),"elasticity"))
    return this->SIMsupel::parse(elem);

  std::string supId;
  if (!utl::getAttribute(elem,"supId",supId))
    return true; // Not a superelement

  // This <elasticity> block defines the FE model of a superelement
  IFEM::cout <<"  Parsing FE model for superelement \""<< supId <<"\"\n";
  SIMgeneric* supEl = new SIMLinEl3D("Linear elastic superelement",false);
  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!supEl->parse(child))
    {
      delete supEl;
      std::cerr <<" *** Failure."<< std::endl;
      return false;
    }

  mySubSim[supId].sim = supEl;
  IFEM::cout <<"  FE model \""<< supId <<"\" loaded."<< std::endl;

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
    if (sub.second.sim)
      ok &= (sub.second.sim->preprocess({},fixDup) &&
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
  if (psol.empty())
    return idBlock;

  VTF* vtf = this->getVTF();
  if (!vtf)
    return idBlock;

  // Recover internal displacements for all superelements
  if (!this->recoverInternalDOFs(psol))
    return -9;

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

int SIMLinElSup::writeGlvS2 (const Vector&, int iStep, int& nBlock,
                             double, int idBlock, int)
{
  VTF* vtf = this->getVTF();
  if (!vtf)
    return idBlock;

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
          return -1;

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
      return -3;
  }

  // Write result block identifications

  if (idBlock <= 20)
    idBlock = 21; // avoid conflict with blocks from the 1D simulator
  bool ok = true;
  for (size_t i = 0; i < sID.size() && !sID[i].empty() && ok; i++)
    ok = vtf->writeSblk(sID[i], problem->getField2Name(i).c_str(),
                        idBlock++, iStep);

  return ok ? idBlock : -5;
}
