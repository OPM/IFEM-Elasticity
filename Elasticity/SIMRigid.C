// $Id$
//==============================================================================
//!
//! \file SIMRigid.C
//!
//! \date Dec 09 2020
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Rigid and nodal coupling handler for elasticity problems.
//!
//==============================================================================

#include "SIMRigid.h"
#include "SIMinput.h"
#include "ASMbase.h"
#include "Utilities.h"
#include "ElementBlock.h"
#include "IFEM.h"
#include "Vec3Oper.h"
#include "tinyxml2.h"
#include <iostream>
#include <fstream>
#include <sstream>


bool SIMRigid::parseRigid (const tinyxml2::XMLElement* elem, SIMinput* mySim)
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


bool SIMRigid::parseCouplings (const tinyxml2::XMLElement* elem)
{
  // Lambda function parsing the nodal couplings for a patch.
  auto&& parseCpl = [](std::istream& is, Couplin& cpl, char delim = '\n')
  {
    std::string cline;
    int slave, master;
    double weight;
    while (std::getline(is,cline,delim))
    {
      std::istringstream ss(cline);
      ss >> slave >> master >> weight;
      while (ss)
      {
        cpl[slave].emplace_back(master,weight);
        ss >> master >> weight;
      }
#if INT_DEBUG > 1
      std::cout <<"Slave ("<< slave <<") coupled to master(s):";
      for (const Master& mst : cpl[slave])
        std::cout <<" ("<< mst.first <<")*"<< mst.second;
      std::cout << std::endl;
#endif
    }
  };

  IFEM::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;

  int patch = 1;
  utl::getAttribute(elem,"patch",patch);
  if (std::string fName; utl::getAttribute(elem,"file",fName))
  {
    if (std::ifstream fs(fName); fs)
      parseCpl(fs,myCouplings[patch]);
    else
    {
      std::cerr <<" SIMRigid::parseCouplings: Unable to open file \""
                << fName  <<"\"."<< std::endl;
      return false;
    }
  }
  else if (elem->FirstChild())
  {
    std::istringstream ss(elem->FirstChild()->Value());
    parseCpl(ss,myCouplings[patch],'\\');
  }
  else
    return false;

  IFEM::cout <<"\tParsed "<< myCouplings[patch].size()
             <<" nodal couplings for Patch "<< patch;
  if (double xtol = 0.0; utl::getAttribute(elem,"xtol",xtol) && xtol > 0.0)
  {
    myGeomTols[patch] = xtol;
    IFEM::cout <<" (tolerance = "<< xtol <<")";
  }

  IFEM::cout << std::endl;
  return true;
}


bool SIMRigid::addRigidMPCs (SIMinput* mySim, int& ngnod) const
{
  for (PropertyVec::const_iterator pit = mySim->begin_prop();
       pit != mySim->end_prop(); ++pit)
    if (pit->pcode == Property::RIGID)
    {
      PntMap::const_iterator mit = myMasters.find(pit->pindx);
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
      if (SIMinput::IdxVec3* mp = mySim->getDiscretePoint(mit->second.item); mp)
        if (ASMbase* pch = mySim->getPatch(pit->patch); pch)
          if (pch->addRigidCpl(pit->lindx,pit->ldim,pit->basis,
                               mp->first,mp->second)) ++ngnod;
    }

  return true;
}


bool SIMRigid::addGeneralCouplings (SIMinput* mySim) const
{
  bool status = true;
  for (const CplMap::value_type& cpl : myCouplings)
    if (ASMbase* pch = mySim->getPatch(cpl.first); pch)
    {
      std::vector<Ipair> connectedNodes;
      for (const Couplin::value_type& mpc : cpl.second)
        if (mpc.second.size() == 1 && mpc.second.front().second == 1.0)
          // Direct coupling of matching nodes
          connectedNodes.emplace_back(mpc.first,mpc.second.front().first);
        else
        {
          // Create linear couplings
          IntVec    masters;
          RealArray weights;
          for (const Master& mst : mpc.second)
          {
            masters.push_back(mst.first);
            weights.push_back(mst.second);
          }
          pch->addNodalCouplings(mpc.first,masters,weights);
        }

      if (TolMap::const_iterator it = myGeomTols.find(cpl.first);
          it != myGeomTols.end())
        status &= pch->selfInterconnect(connectedNodes,it->second);
      else
        status &= pch->selfInterconnect(connectedNodes);
    }

  return status;
}


ElementBlock* SIMRigid::rigidGeometry (SIMinput* mySim) const
{
  if (myMasters.empty()) return nullptr; // no rigid couplings

  size_t inod = 0;
  std::map<int,size_t> imap;
  ElementBlock* rgd = new ElementBlock(2);
  rgd->unStructResize(0,myMasters.size());
  for (const PntMap::value_type& master : myMasters)
  {
    if (SIMinput::IdxVec3* mp = mySim->getDiscretePoint(master.second.item); mp)
      rgd->setCoor(inod,mp->second);
    imap[master.first] = inod++;
  }

  int slvThick = mySim->opt.discretization == ASM::SplineC1 ? 2 : 1;

  for (PropertyVec::const_iterator pit = mySim->begin_prop();
       pit != mySim->end_prop(); ++pit)
    if (pit->pcode == Property::RIGID)
      if (PntMap::const_iterator mit = myMasters.find(pit->pindx);
          mit != myMasters.end() && mit->second.patch == 0)
        if (ASMbase* pch = mySim->getPatch(pit->patch); pch &&
            mySim->getDiscretePoint(mit->second.item))
        {
          IntVec nodes;
          pch->getBoundaryNodes(pit->lindx,nodes,pit->basis,slvThick,0,true);
          for (int node : nodes)
            rgd->addLine(imap[mit->first],pch->getCoord(node));
        }

  return rgd;
}
