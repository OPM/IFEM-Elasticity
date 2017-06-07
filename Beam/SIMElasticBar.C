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


SIMElasticBar::SIMElasticBar (unsigned char n) : SIM1D(3)
{
  nsd = 3;
  nsv = n;
  printed = false;
}


SIMElasticBar::~SIMElasticBar ()
{
  for (LoadMap::iterator it = myLoads.begin(); it != myLoads.end(); ++it)
    delete it->second;
}


void SIMElasticBar::printProblem() const
{
  if (printed) return;

  this->SIM1D::printProblem();
  printed = true; // Avoid printing problem definition more than once
}


ElasticBar* SIMElasticBar::getBarIntegrand (const std::string& type)
{
  if (type == "cable")
    myProblem = new ElasticCable(nsv);
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

  if (!bar && !beam) return false;

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

    else if (!strcasecmp(child->Value(),"nodeload"))
    {
      int patch = 1, node = 0, dof = 0;
      std::string type("constant");
      utl::getAttribute(child,"patch",patch);
      utl::getAttribute(child,"node",node);
      utl::getAttribute(child,"dof",dof);
      utl::getAttribute(child,"type",type);

      Vec3 Xnod;
      double xi;
      if (utl::getAttribute(child,"u",xi) && xi >= 0.0 && xi <= 1.0)
      {
        for (PatchVec::iterator it = myModel.begin(); it != myModel.end(); ++it)
          (*it)->setNoFields(nf);
        if (this->createFEMmodel())
          node = this->evalPoint(&xi,Xnod,&xi,patch,true);
      }

      if (child->FirstChild() && node > 0 && dof > 0 && dof <= nf)
      {
        ScalarFunc* f = nullptr;
        IFEM::cout <<"\tNode "<< node <<" dof "<< dof <<" Load: ";
        if (type == "constant")
        {
          f = new ConstantFunc(atof(child->FirstChild()->Value()));
          IFEM::cout << (*f)(0.0) << std::endl;
        }
        else
          f = utl::parseTimeFunc(child->FirstChild()->Value(),type);

        myLoads[std::make_pair(node,dof)] = f;
      }
    }
    else
      myProblem->parse(child);
  }

  return true;
}


void SIMElasticBar::preprocessA ()
{
  this->printProblem();

  for (PatchVec::iterator it = myModel.begin(); it != myModel.end(); ++it)
    (*it)->setNoFields(nf);
}


bool SIMElasticBar::renumberNodes (const std::map<int,int>& nodeMap)
{
  bool ok = this->SIM1D::renumberNodes(nodeMap);

  LoadMap newLoads;
  for (LoadMap::iterator it = myLoads.begin(); it != myLoads.end(); ++it)
  {
    int node = it->first.first;
    if (utl::renumber(node,nodeMap,true))
      newLoads[std::make_pair(node,it->first.second)] = it->second;
    else
      ok = false;
  }

  std::swap(myLoads,newLoads);
  return ok;
}


bool SIMElasticBar::assembleDiscreteTerms (const IntegrandBase* itg,
                                           const TimeDomain& time)
{
  if (itg != myProblem)
    return true;

  SystemVector* R = myEqSys->getVector(2); // External load vector
  double scale = 1.0;

  if (!R || itg->getIntegrationPrm(4) != 1.0)
  {
    R = myEqSys->getVector(0); // System right-hand-side vector
    if (itg->getIntegrationPrm(3) <= 0.0) // HHT is used
      scale = itg->getIntegrationPrm(2) + 1.0; // alphaH + 1.0
  }
  if (!R) return true; // Silently ignore, if no right-hand-side vector

  // Assemble external nodal point loads at current time step
  for (LoadMap::const_iterator it = myLoads.begin(); it != myLoads.end(); ++it)
    if (!mySam->assembleSystem(*R,(*it->second)(time.t)*scale,it->first))
      return false;

  return true;
}


Tensor SIMElasticBar::getNodeRotation (int inod) const
{
  for (PatchVec::const_iterator it = myModel.begin(); it != myModel.end(); ++it) {
    size_t node;
    if ((node = (*it)->getNodeIndex(inod,true)))
      return static_cast<const ASMs1D*>(*it)->getRotation(node);
  }

  return Tensor(nsd,true);
}
