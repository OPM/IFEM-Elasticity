// $Id$
//==============================================================================
//!
//! \file SIMRigid.C
//!
//! \date Dec 09 2020
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Rigid Coupling handler for elasticity problems.
//!
//==============================================================================

#include "SIMRigid.h"
#include "SIMinput.h"
#include "ASMbase.h"
#include "Vec3Oper.h"
#include "Utilities.h"
#include "ElementBlock.h"
#include "IFEM.h"
#include "tinyxml.h"


bool SIMRigid::parseRigid (const TiXmlElement* elem, SIMinput* mySim)
{
  IFEM::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;

  std::string master, slave;
  utl::getAttribute(elem,"set",slave);
  utl::getAttribute(elem,"slave",slave);
  utl::getAttribute(elem,"master",master);

  int islave = mySim->getUniquePropertyCode(slave);
  if (islave == 0) return false;

  const TopEntity& topEnt = mySim->getEntity(master);
  if (topEnt.empty())
  {
    std::cerr <<" *** SIMRigid::parse: Undefined/empty topology set \""
              << master <<"\"."<< std::endl;
    return false;
  }

  TopEntity::const_iterator titem = topEnt.begin();
  if (topEnt.size() != 1 || titem->idim > 0)
  {
    std::cerr <<" *** SIMRigid::parse: Invalid topology set \""
              << master <<"\". It should contain a single point."<< std::endl;
    return false;
  }

  IFEM::cout <<"\tSlave code "<< islave <<" ("<< master <<"):";
  if (titem->patch > 0)
    IFEM::cout <<" Patch/item index "<< titem->patch <<", "<< titem->item;
  else
    IFEM::cout <<" Master point index "<< titem->item
               <<" ("<< mySim->getDiscretePoint(titem->item)->second <<")";
  IFEM::cout << std::endl;

  mySim->setPropertyType(islave,Property::RIGID);
  myMasters[islave] = *titem;
  return true;
}


bool SIMRigid::addRigidMPCs (SIMinput* mySim, int& ngnod) const
{
  PropertyVec::const_iterator pit;
  for (pit = mySim->begin_prop(); pit != mySim->end_prop(); ++pit)
    if (pit->pcode == Property::RIGID)
    {
      std::map<int,TopItem>::const_iterator mit = myMasters.find(pit->pindx);
      if (mit == myMasters.end())
      {
        std::cerr <<" *** SIMRigid::addRigidMPCs: Local error"
                  <<", unknown rigid code "<< pit->pindx << std::endl;
        return false;
      }
      else if (mit->second.patch > 0)
      {
        IFEM::cout <<" *** SIMRigid::addRigidMPCs:"
                   <<" Sorry, patch vertex as master not implemented yet."
                   << std::endl;
        return false;
      }

      // Find master point and the affected patch,
      // and create a rigid coupling between them
      SIMinput::IdxVec3* mst = mySim->getDiscretePoint(mit->second.item);
      ASMbase* pch = mySim->getPatch(pit->patch);
      if (mst && pch)
        if (pch->addRigidCpl(pit->lindx,pit->ldim,pit->basis,
                             mst->first,mst->second)) ++ngnod;
    }

  return true;
}


ElementBlock* SIMRigid::rigidGeometry (SIMinput* mySim) const
{
  if (myMasters.empty()) return nullptr; // no rigid couplings

  size_t inod = 0;
  std::map<int,size_t> imap;
  ElementBlock* rgd = new ElementBlock(2);
  rgd->unStructResize(0,myMasters.size());
  for (const std::pair<int,TopItem>& master : myMasters)
  {
    SIMinput::IdxVec3* mst = mySim->getDiscretePoint(master.second.item);
    if (mst) rgd->setCoor(inod,mst->second);
    imap[master.first] = inod++;
  }

  PropertyVec::const_iterator pit;
  std::map<int,TopItem>::const_iterator mit;
  for (pit = mySim->begin_prop(); pit != mySim->end_prop(); ++pit)
    if (pit->pcode == Property::RIGID)
      if ((mit = myMasters.find(pit->pindx)) != myMasters.end())
        if (mit->second.patch == 0)
        {
          SIMinput::IdxVec3* mst = mySim->getDiscretePoint(mit->second.item);
          ASMbase* pch = mySim->getPatch(pit->patch);
          if (mst && pch)
          {
            IntVec nodes;
            pch->getBoundaryNodes(pit->lindx,nodes,pit->basis,1,0,true);
            for (int node : nodes)
              rgd->addLine(imap[mit->first],pch->getCoord(node));
          }
        }

  return rgd;
}
