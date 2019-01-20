// $Id$
//==============================================================================
//!
//! \file SIMElasticBar.C
//!
//! \date Aug 11 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based FE analysis of elastic bars & beams.
//!
//==============================================================================

#include "SIMElasticBar.h"
#include "ElasticCable.h"
#include "ElasticBeam.h"
#include "AlgEqSystem.h"
#include "ASMs1D.h"
#include "SAM.h"
#include "AnaSol.h"
#include "Functions.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include "tinyxml.h"


SIMElasticBar::SIMElasticBar (unsigned char n) : SIMElastic1D(3)
{
  nsd = 3;
  nsv = n;
  printed = false;
}


SIMElasticBar::~SIMElasticBar ()
{
  for (PointLoad& load : myLoads)
    delete load.p;
}


void SIMElasticBar::printProblem() const
{
  if (!printed)
    this->SIM1D::printProblem();

  printed = true; // Avoid printing problem definition more than once
}


ElasticBar* SIMElasticBar::getBarIntegrand (const std::string& type)
{
  if (type == "cable")
    myProblem = new ElasticCable(nsd,nsv);
  else
    myProblem = new ElasticBar(toupper(type[0]),nsd);

  // Always one integration point per element
  if (type != "cable")
    opt.nGauss[0] = opt.nGauss[1] = 1;

  return dynamic_cast<ElasticBar*>(myProblem);
}


ElasticBeam* SIMElasticBar::getBeamIntegrand (const std::string&)
{
  myProblem = new ElasticBeam(nsv);

  // Always one integration point per element
  opt.nGauss[0] = opt.nGauss[1] = 1;

  return dynamic_cast<ElasticBeam*>(myProblem);
}


bool SIMElasticBar::parse (const TiXmlElement* elem)
{
  bool isCable = !strcasecmp(elem->Value(),"cable");
  if (isCable || !strcasecmp(elem->Value(),"bar"))
    nf = 3;
  else if (!strcasecmp(elem->Value(),"beam"))
    nf = 6;
  else if (!strcasecmp(elem->Value(),"anasol"))
  {
    std::string type("expression");
    utl::getAttribute(elem,"type",type,true);
    if (type == "expression")
    {
      IFEM::cout <<"\tAnalytical solution: Expression"<< std::endl;
      if (!mySol) mySol = new AnaSol(elem,false);
    }
    else
      std::cerr <<"  ** SIMElasticBar::parse: Invalid analytical solution "
                << type <<" (ignored)"<< std::endl;
    return true;
  }
  else
    return this->SIM1D::parse(elem);

  ElasticBar*  bar  = nullptr;
  ElasticBeam* beam = nullptr;
  if (myProblem)
    switch (nf) {
    case 3:
      bar = dynamic_cast<ElasticBar*>(myProblem);
      break;
    case 6:
      beam = dynamic_cast<ElasticBeam*>(myProblem);
      break;
    }
  else
  {
    std::string type(isCable ? "cable" : "N");
    utl::getAttribute(elem,"type",type,true);
    switch (nf) {
    case 3:
      bar = this->getBarIntegrand(type);
      break;
    case 6:
      beam = this->getBeamIntegrand(type);
      break;
    }
  }

  if (!bar && !beam)
    return false;

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
  {
    IFEM::cout <<"  Parsing <"<< child->Value() <<">"<< std::endl;
    if (!strcasecmp(child->Value(),"gravity"))
    {
      Vec3 g;
      utl::getAttribute(child,"x",g.x);
      utl::getAttribute(child,"y",g.y);
      utl::getAttribute(child,"z",g.z);
      IFEM::cout <<"\tGravitation vector: "<< g << std::endl;
      if (bar)
        bar->setGravity(g);
      else
        beam->setGravity(g);
    }
    else if (!strcasecmp(child->Value(),"material"))
    {
      double E = 2.1e11, G = 8.1e10, nu = 0.0, rho = 7.85e3;
      utl::getAttribute(child,"E",E);
      if (!utl::getAttribute(child,"G",G) && utl::getAttribute(child,"nu",nu))
        if (nu >= 0.0 && nu < 0.5)
          G = E/(nu+nu+2.0); // Derive G-modulus from E and nu
      utl::getAttribute(child,"rho",rho);
      if (bar)
      {
        double A = 0.1, I = 0.0;
        ElasticBeam::parsePipe(child,A,I); // "D" or "R" and optionally "t"
        double EA = 0.0, EI = 0.0;
        if (!utl::getAttribute(child,"EA",EA))
        {
          utl::getAttribute(child,"A",A);
          EA = E*A;
        }
        bar->setStiffness(EA);
        if (isCable)
        {
          if (!utl::getAttribute(child,"EI",EI))
          {
            utl::getAttribute(child,"I",I);
            EI = E*I;
          }
          static_cast<ElasticCable*>(bar)->setBendingStiffness(EI);
        }
        bar->setMass(rho);
        IFEM::cout <<"\tAxial stiffness = "<< EA;
        if (isCable && EI > 0.0)
          IFEM::cout <<"\n\tBending stiffness = "<< EI;
      }
      else
      {
        beam->setStiffness(E,G);
        beam->setMass(rho);
        IFEM::cout <<"\tStiffness moduli = "<< E <<" "<< G;
      }
      IFEM::cout << (isCable ? "\n\tM" : ", m")
                  <<"ass density = "<< rho << std::endl;
    }
    else if (beam && !strcasecmp(child->Value(),"properties"))
      beam->parseBeamProperties(child);

    else if (beam && !strcasecmp(child->Value(),"twist"))
      this->parseTwist(child);

    else if (beam && !strcasecmp(child->Value(),"lineload"))
      beam->parseBeamLoad(child);

    else if (beam && !strcasecmp(child->Value(),"cplload"))
      beam->parseBeamLoad(child);

    else if (!strcasecmp(child->Value(),"pointload") && child->FirstChild())
    {
      PointLoad load(1);
      utl::getAttribute(child,"dof",load.ldof);
      utl::getAttribute(child,"u",load.xi);
      if (load.ldof > 0 && load.ldof <= nf && load.xi >= 0.0 && load.xi <= 1.0)
      {
        // Negate the local DOF flag to signal that element loads are allowed
        bool allowElementPointLoad = false;
        utl::getAttribute(child,"onElement",allowElementPointLoad);
        if (allowElementPointLoad) load.ldof = -load.ldof;

        if (utl::getAttribute(child,"patch",load.inod))
          IFEM::cout <<"\tPoint: P"<< load.inod;
        else
          IFEM::cout <<"\tPoint:";
        IFEM::cout <<" xi = "<< load.xi <<" dof = "<< load.ldof <<" Load: ";

        std::string type("constant");
        utl::getAttribute(child,"type",type);
        if (type == "constant")
        {
          load.p = new ConstantFunc(atof(child->FirstChild()->Value()));
          IFEM::cout << (*load.p)(0.0) << std::endl;
        }
        else
          load.p = utl::parseTimeFunc(child->FirstChild()->Value(),type);

        myLoads.push_back(load);
      }
    }

    else if (!strcasecmp(child->Value(),"nodeload") && child->FirstChild())
    {
      PointLoad load;
      utl::getAttribute(child,"node",load.inod);
      utl::getAttribute(child,"dof",load.ldof);

      if (load.inod > 0 && load.ldof > 0 && load.ldof <= nf)
      {
        std::string type("constant");
        utl::getAttribute(child,"type",type);

        IFEM::cout <<"\tNode "<< load.inod <<" dof "<< load.ldof <<" Load: ";
        if (type == "constant")
        {
          load.p = new ConstantFunc(atof(child->FirstChild()->Value()));
          IFEM::cout << (*load.p)(0.0) << std::endl;
        }
        else
          load.p = utl::parseTimeFunc(child->FirstChild()->Value(),type);

        myLoads.push_back(load);
      }
    }
    else
      myProblem->parse(child);
  }

  return true;
}


void SIMElasticBar::preprocessA ()
{
  if (nf == 3 && nsd < 3)
    nf = nsd; // 2D bar/cable, two DOFs per node

  this->printProblem();
  for (ASMbase* pch : myModel)
    pch->setNoFields(nf);
}


bool SIMElasticBar::preprocessB ()
{
  // Preprocess the nodal point loads, if any
  if (myLoads.empty())
    return true;

  int ipt = 0;
  bool ok = true;
  for (PointLoad& pl : myLoads)
    if (pl.xi >= 0.0 && pl.inod > 0)
    {
      int foundPoint = this->findLoadPoint(++ipt,pl.inod,pl.xi,pl.ldof < 0);
      if (foundPoint > 0)
        pl.inod = foundPoint;
      else if (foundPoint == 0)
        ok = false;
    }

  if (ipt > 0)
    IFEM::cout <<"\n"<< std::endl;

  return ok;
}


bool SIMElasticBar::renumberNodes (const std::map<int,int>& nodeMap)
{
  bool ok = this->SIM1D::renumberNodes(nodeMap);

  for (PointLoad& load : myLoads)
    if (load.xi < 0.0 && load.inod > 0)
      ok &= utl::renumber(load.inod,nodeMap,true);

  return ok;
}


bool SIMElasticBar::assembleDiscreteTerms (const IntegrandBase* itg,
                                           const TimeDomain& time)
{
  bool ok = true;
  if (itg != myProblem)
    return ok;

  SystemVector* R = myEqSys->getVector(2); // External load vector
  double scl = 1.0;

  SIM::SolutionMode mode = itg->getMode();
  if (!myLoads.empty() && mode == SIM::ARCLEN)
    this->setMode(SIM::RHS_ONLY);

  if (!R || itg->getIntegrationPrm(4) != 1.0)
  {
    R = myEqSys->getVector(0); // System right-hand-side vector
    if (itg->getIntegrationPrm(3) <= 0.0) // HHT is used
      scl = itg->getIntegrationPrm(2) + 1.0; // alphaH + 1.0
  }

  if (R)
  {
    // Assemble external nodal point loads at current time step
    for (const PointLoad& load : myLoads)
      if (load.ldof > 0)
        ok &= mySam->assembleSystem(*R,(*load.p)(time.t)*scl,
                                    std::make_pair(load.inod,load.ldof));
      else // This is an element point load
        ok &= this->assemblePoint(load.inod,load.xi,(*load.p)(time.t),
                                  -load.ldof);
  }

  if (mode == SIM::ARCLEN)
    R = myEqSys->getVector(1); // External load gradient for arc-length driver
  else
    R = nullptr;

  if (R)
  {
    // Assemble external nodal point load gradient at current time step
    this->setMode(mode);
    static_cast<ElasticBase*>(myProblem)->setLoadGradientMode();
    for (const PointLoad& load : myLoads)
      if (load.ldof > 0)
        ok &= mySam->assembleSystem(*R,load.p->deriv(time.t),
                                    std::make_pair(load.inod,load.ldof));
      else // This is an element point load
        ok &= this->assemblePoint(load.inod,load.xi,load.p->deriv(time.t),
                                  -load.ldof);
  }

  return ok;
}


Tensor SIMElasticBar::getNodeRotation (int inod) const
{
  for (const ASMbase* pch : myModel)
  {
    size_t node = pch->getNodeIndex(inod,true);
    if (node > 0)
      return static_cast<const ASMs1D*>(pch)->getRotation(node);
  }

  return Tensor(nsd,true);
}


bool SIMElasticBar::getExtLoad (RealArray& extloa, const TimeDomain& time) const
{
  extloa.resize(nf,0.0);
  for (size_t i = 0; i < nf; i++)
    extloa[i] = this->extractScalar(i);

  for (const PointLoad& load : myLoads)
    if (load.ldof > 0 && load.ldof <= nf)
      extloa[load.ldof-1] += (*load.p)(time.t);

  return true;
}
