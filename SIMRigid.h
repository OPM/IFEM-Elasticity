// $Id$
//==============================================================================
//!
//! \file SIMRigid.h
//!
//! \date Dec 09 2020
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Rigid Coupling handler for elasticity problems.
//!
//=============================================================================

#ifndef _SIM_RIGID_H
#define _SIM_RIGID_H

#include "TopologySet.h"

class SIMinput;
class ElementBlock;
class TiXmlElement;


/*!
  \brief Rigid coupling handler for elasticity problems.
*/

class SIMRigid
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  SIMRigid() {}
  //! \brief Empty destructor.
  virtual ~SIMRigid() {}

  //! \brief Parses a rigid coupling definition from an XML element.
  bool parseRigid(const TiXmlElement* elem, SIMinput* mySim);

  //! \brief Creates multi-point constraint equations for the rigid couplings.
  bool addRigidMPCs(SIMinput* mySim, int& ngnod) const;

  //! \brief Creates an element block visualizing the rigid couplings.
  ElementBlock* rigidGeometry(SIMinput* mySim) const;

private:
  std::map<int,TopItem> myMasters; //!< Discrete master points
};

#endif
