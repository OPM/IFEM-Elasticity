// $Id$
//==============================================================================
//!
//! \file SIMRigid.h
//!
//! \date Dec 09 2020
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Rigid and nodal coupling handler for elasticity problems.
//!
//=============================================================================

#ifndef _SIM_RIGID_H
#define _SIM_RIGID_H

#include "TopologySet.h"
#include <vector>

class SIMinput;
class ElementBlock;
namespace tinyxml2 { class XMLElement; }


/*!
  \brief Rigid and nodal coupling handler for elasticity problems.
*/

class SIMRigid
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  SIMRigid() {}
  //! \brief Empty destructor.
  virtual ~SIMRigid() {}

  //! \brief Parses a rigid coupling definition from an XML element.
  bool parseRigid(const tinyxml2::XMLElement* elem, SIMinput* mySim);
  //! \brief Parses general nodal couplings from an XML element.
  bool parseCouplings(const tinyxml2::XMLElement* elem);

  //! \brief Creates multi-point constraint equations for the rigid couplings.
  bool addRigidMPCs(SIMinput* mySim, int& ngnod) const;
  //! \brief Creates multi-point constraint equations for the general couplings.
  bool addGeneralCouplings(SIMinput* mySim) const;

  //! \brief Creates an element block visualizing the rigid couplings.
  ElementBlock* rigidGeometry(SIMinput* mySim) const;

private:
  using PntMap  = std::map<int,TopItem>; //!< Discrete master point definition
  using Master  = std::pair<int,double>; //!< Independent node with weight
  using Masters = std::vector<Master>;   //!< Independent node list
  using Couplin = std::map<int,Masters>; //!< General nodal coupling
  using CplMap  = std::map<int,Couplin>; //!< Couplings for patches
  using TolMap  = std::map<int,double>;  //!< Geometry tolerance for patches

  PntMap myMasters;   //!< Discrete master point container
  CplMap myCouplings; //!< Nodal coupling container
  TolMap myGeomTols;  //!< Geometry tolerance container
};

#endif
