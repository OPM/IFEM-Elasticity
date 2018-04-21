// $Id$
//==============================================================================
//!
//! \file SIMLinElKL.C
//!
//! \date Sep 16 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based FE analysis of Kirchhoff-Love shells.
//!
//==============================================================================

#include "SIMLinElKL.h"
#include "KirchhoffLovePlate.h"
#include "KirchhoffLoveShell.h"
#include "AnalyticSolutions.h"
#include "AnaSol.h"
#include "tinyxml.h"


KirchhoffLove* SIMLinElKL::getProblem (int version)
{
  KirchhoffLove* klp = dynamic_cast<KirchhoffLove*>(myProblem);
  if (!myProblem)
  {
    if (nsd == 3)
      myProblem = klp = new KirchhoffLoveShell();
    else
      myProblem = klp = new KirchhoffLovePlate(2,version);
  }

  return klp;
}


bool SIMLinElKL::parse (char* keyWord, std::istream& is)
{
  if (!strncasecmp(keyWord,"ANASOL",6))
  {
    char* cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"NAVIERPLATE",7))
    {
      double a  = atof(strtok(NULL," "));
      double b  = atof(strtok(NULL," "));
      double t  = atof(strtok(NULL," "));
      double E  = atof(strtok(NULL," "));
      double nu = atof(strtok(NULL," "));
      double pz = atof(strtok(NULL," "));
      IFEM::cout <<"\nAnalytic solution: NavierPlate a="<< a <<" b="<< b
                 <<" t="<< t <<" E="<< E <<" nu="<< nu <<" pz="<< pz;
      if ((cline = strtok(NULL," ")))
      {
	double xi  = atof(cline);
	double eta = atof(strtok(NULL," "));
	IFEM::cout <<" xi="<< xi <<" eta="<< eta;
	if ((cline = strtok(NULL," ")))
	{
	  double c = atof(cline);
	  double d = atof(strtok(NULL," "));
	  IFEM::cout <<" c="<< c <<" d="<< d;
	  if (!mySol)
	    mySol = new NavierPlate(a,b,t,E,nu,pz,xi,eta,c,d);
	}
	else if (!mySol)
	  mySol = new NavierPlate(a,b,t,E,nu,pz,xi,eta);
      }
      else if (!mySol)
	mySol = new NavierPlate(a,b,t,E,nu,pz);
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      IFEM::cout <<"\nAnalytical solution: Expression"<< std::endl;
      int lines = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      if (!mySol)
        mySol = new AnaSol(is,lines,false);
    }
    else
    {
      std::cerr <<"  ** SIMLinElKL::parse: Unknown analytical solution "
                << cline <<" (ignored)"<< std::endl;
      return true;
    }
  }

  else
    return this->SIMKLShell::parse(keyWord,is);

  return true;
}


bool SIMLinElKL::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"kirchhofflove"))
    return this->SIM2D::parse(elem);
  else if (!this->SIMKLShell::parse(elem))
    return false;

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(),"anasol"))
    {
      std::string type;
      utl::getAttribute(child,"type",type,true);
      if (type == "navierplate") {
        double a = 0.0, b = 0.0, c = 0.0, d = 0.0, t = 0.0;
        double E = 10000.0, nu = 0.3, pz = 1.0, xi = 0.0, eta = 0.0;
        int max_mn = 100;
        utl::getAttribute(child,"a",a);
        utl::getAttribute(child,"b",b);
        utl::getAttribute(child,"c",c);
        utl::getAttribute(child,"d",d);
        utl::getAttribute(child,"t",t);
        utl::getAttribute(child,"E",E);
        utl::getAttribute(child,"nu",nu);
        utl::getAttribute(child,"pz",pz);
        utl::getAttribute(child,"xi",xi);
        utl::getAttribute(child,"eta",eta);
        utl::getAttribute(child,"nTerm",max_mn);
        IFEM::cout <<"\tAnalytic solution: NavierPlate a="<< a <<" b="<< b
                   <<" t="<< t <<" E="<< E <<" nu="<< nu <<" pz="<< pz;
        if (xi != 0.0 && eta != 0.0) {
          IFEM::cout <<" xi="<< xi <<" eta="<< eta;
          if (c != 0.0 && d != 0.0) {
            IFEM::cout <<" c="<< c <<" d="<< d;
            if (!mySol)
              mySol = new NavierPlate(a,b,t,E,nu,pz,xi,eta,c,d,max_mn);
          }
          else if (!mySol)
            mySol = new NavierPlate(a,b,t,E,nu,pz,xi,eta,0.0,0.0,max_mn);
        }
        else if (!mySol)
          mySol = new NavierPlate(a,b,t,E,nu,pz,max_mn);
      }
      else if (type == "circularplate") {
        double R = 0.0, t = 0.0, E = 10000.0, nu = 0.3, P = 1000000.0;
        utl::getAttribute(child,"R",R);
        utl::getAttribute(child,"t",t);
        utl::getAttribute(child,"E",E);
        utl::getAttribute(child,"nu",nu);
        utl::getAttribute(child,"P",P);
        IFEM::cout <<"\tAnalytic solution: CircularPlate R="<< R <<" t="<< t
                   <<" E="<< E <<" nu="<< nu <<" P="<< P << std::endl;
        if (!mySol)
          mySol = new CircularPlate(R,t,E,nu,P);
      }
      else if (type == "expression") {
        IFEM::cout <<"\tAnalytical solution: Expression"<< std::endl;
        if (!mySol)
          mySol = new AnaSol(child);
      }
      else
        std::cerr <<"  ** SIMLinElKL::parse: Unknown analytical solution "
                  << type <<" (ignored)"<< std::endl;
    }

  return true;
}


void SIMLinElKL::preprocessA ()
{
  this->getProblem();
  this->printProblem();

  ThinPlateSol* plSol = dynamic_cast<ThinPlateSol*>(mySol);
  if (!plSol) return;

  // Define analytical boundary condition fields (for rotations)
  for (Property& prop : myProps)
    if (prop.pcode == Property::DIRICHLET_ANASOL)
    {
      if (abs(prop.pindx) >= 200)
      {
        if (aCode[2] == abs(prop.pindx))
          prop.pcode = Property::DIRICHLET_INHOM;
        else if (aCode[2] == 0 && plSol->thetaY())
        {
          aCode[2] = abs(prop.pindx);
          myScalars[aCode[2]] = plSol->thetaY();
          prop.pcode = Property::DIRICHLET_INHOM;
        }
        else
          prop.pcode = Property::UNDEFINED;
      }
      else if (abs(prop.pindx) >= 100)
      {
        if (aCode[1] == abs(prop.pindx))
          prop.pcode = Property::DIRICHLET_INHOM;
        else if (aCode[1] == 0 && plSol->thetaX())
        {
          aCode[1] = abs(prop.pindx);
          myScalars[aCode[1]] = plSol->thetaX();
          prop.pcode = Property::DIRICHLET_INHOM;
        }
        else
          prop.pcode = Property::UNDEFINED;
      }
      else if (abs(prop.pindx) > 0)
      {
        if (aCode[0] == abs(prop.pindx))
          prop.pcode = Property::DIRICHLET_INHOM;
        else if (aCode[0] == 0 && mySol->getScalarSol())
        {
          aCode[0] = abs(prop.pindx);
          myScalars[aCode[0]] = plSol->getScalarSol();
          prop.pcode = Property::DIRICHLET_INHOM;
        }
        else
          prop.pcode = Property::UNDEFINED;
      }
    }
}


bool SIMLinElKL::haveAnaSol () const
{
  if (!mySol) return false;

  KirchhoffLovePlate* klp = dynamic_cast<KirchhoffLovePlate*>(myProblem);
  if (!klp) return false;

  if (klp->getVersion() == 1 && mySol->getStressSol())
    return true;
  else if (klp->getVersion() > 1 && mySol->getScalarSecSol())
    return true;
  else
    return false;
}
