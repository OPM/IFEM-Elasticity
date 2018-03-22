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
#include "LinIsotropic.h"
#include "AnalyticSolutions.h"
#include "AlgEqSystem.h"
#include "ASMbase.h"
#include "SAM.h"
#include "Functions.h"
#include "Utilities.h"
#include "Property.h"
#include "AnaSol.h"
#include "Vec3Oper.h"
#include "tinyxml.h"


SIMLinElKL::SIMLinElKL (bool shell)
{
  if (shell) nsd = 3;

  nf[0] = shell ? 3 : 1;
  aCode[0] = aCode[1] = aCode[2] = 0;
}


SIMLinElKL::~SIMLinElKL ()
{
  // To prevent the SIMbase destructor try to delete already deleted functions
  for (int i = 0; i < 3; i++)
    if (aCode[i] > 0) myScalars.erase(aCode[i]);
}


Elasticity* SIMLinElKL::getIntegrand ()
{
  return dynamic_cast<Elasticity*>(myProblem);
}


bool SIMLinElKL::parse (char* keyWord, std::istream& is)
{
  KirchhoffLove* klp = dynamic_cast<KirchhoffLove*>(myProblem);
  if (!klp)
  {
    if (nsd == 3)
      myProblem = klp = new KirchhoffLoveShell();
    else
      myProblem = klp = new KirchhoffLovePlate();
  }

  char* cline;
  if (!strncasecmp(keyWord,"GRAVITY",7))
  {
    double g = atof(strtok(keyWord+7," "));
    klp->setGravity(g);
    IFEM::cout <<"\nGravitation constant: "<< g << std::endl;
  }

  else if (!strncasecmp(keyWord,"ISOTROPIC",9))
  {
    int nmat = atoi(keyWord+10);
    IFEM::cout <<"\nNumber of isotropic materials: "<< nmat << std::endl;

    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      if (code > 0)
        this->setPropertyType(code,Property::MATERIAL,mVec.size());

      double E   = atof(strtok(NULL," "));
      double nu  = atof(strtok(NULL," "));
      double rho = (cline = strtok(NULL, " ")) ? atof(cline) : 0.0;
      double thk = (cline = strtok(NULL, " ")) ? atof(cline) : 0.0;
      mVec.push_back(new LinIsotropic(E,nu,rho,true));
      tVec.push_back(thk);
      IFEM::cout <<"\tMaterial code "<< code <<": "
                 << E <<" "<< nu <<" "<< rho <<" "<< thk << std::endl;
    }

    if (!mVec.empty())
      klp->setMaterial(mVec.front());
    if (!tVec.empty() && tVec.front() != 0.0)
      klp->setThickness(tVec.front());
  }

  else if (!strncasecmp(keyWord,"POINTLOAD",9))
  {
    int nload = atoi(keyWord+9);
    IFEM::cout <<"\nNumber of point loads: "<< nload;

    myLoads.resize(nload);
    for (int i = 0; i < nload && (cline = utl::readLine(is)); i++)
    {
      myLoads[i].patch = atoi(strtok(cline," "));
      myLoads[i].xi[0] = atof(strtok(NULL," "));
      myLoads[i].xi[1] = atof(strtok(NULL," "));
      myLoads[i].pload = atof(strtok(NULL," "));
      IFEM::cout <<"\n\tPoint "<< i+1 <<": P"<< myLoads[i].patch
                 <<" xi = "<< myLoads[i].xi[0] <<" "<< myLoads[i].xi[1];
      if (nsd == 3)
      {
        myLoads[i].ldof = (cline = strtok(NULL, " ")) ? atoi(cline) : 1;
        IFEM::cout <<" direction = "<< (int)myLoads[i].ldof;
      }
      IFEM::cout <<" load = "<< myLoads[i].pload;
    }
  }

  else if (!strncasecmp(keyWord,"PRESSURE",8))
  {
    int npres = atoi(keyWord+8);
    IFEM::cout <<"\nNumber of pressures: "<< npres << std::endl;

    for (int i = 0; i < npres && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      double p = atof(strtok(NULL," "));
      IFEM::cout <<"\tPressure code "<< code <<": ";
      cline = strtok(NULL," ");
      myScalars[code] = const_cast<RealFunc*>(utl::parseRealFunc(cline,p));
      IFEM::cout << std::endl;
      if (code > 0)
        this->setPropertyType(code,Property::BODYLOAD);
    }
  }

  else if (!strncasecmp(keyWord,"ANASOL",6))
  {
    cline = strtok(keyWord+6," ");
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
    return this->SIM2D::parse(keyWord,is);

  return true;
}


bool SIMLinElKL::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"kirchhofflove"))
    return this->SIM2D::parse(elem);

  KirchhoffLove* klp = dynamic_cast<KirchhoffLove*>(myProblem);
  if (!klp)
  {
    int version = 1;
    utl::getAttribute(elem,"version",version);
    if (nsd == 3)
      myProblem = klp = new KirchhoffLoveShell();
    else
      myProblem = klp = new KirchhoffLovePlate(2,version);
  }

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())

    if (!strcasecmp(child->Value(),"gravity")) {
      double g = 0.0;
      utl::getAttribute(child,"g",g);
      klp->setGravity(g);
      IFEM::cout <<"\tGravitation constant: "<< g << std::endl;
    }

    else if (!strcasecmp(child->Value(),"isotropic")) {
      int code = this->parseMaterialSet(child,mVec.size());

      double E = 1000.0, nu = 0.3, rho = 1.0, thk = 0.1;
      utl::getAttribute(child,"E",E);
      utl::getAttribute(child,"nu",nu);
      utl::getAttribute(child,"rho",rho);
      utl::getAttribute(child,"thickness",thk);

      mVec.push_back(new LinIsotropic(E,nu,rho,true));
      tVec.push_back(thk);
      IFEM::cout <<"\tMaterial code "<< code <<": "
                 << E <<" "<< nu <<" "<< rho <<" "<< thk << std::endl;
      klp->setMaterial(mVec.front());
      if (tVec.front() != 0.0)
        klp->setThickness(tVec.front());
    }

    else if (!strcasecmp(child->Value(),"pointload")) {
      PointLoad load;
      utl::getAttribute(child,"patch",load.patch);
      if (child->FirstChild()) {
        load.pload = atof(child->FirstChild()->Value());
        utl::getAttribute(child,"xi",load.xi[0]);
        utl::getAttribute(child,"eta",load.xi[1]);
        IFEM::cout <<"\tPoint: P"<< load.patch
                   <<" xi = "<< load.xi[0] <<" "<< load.xi[1];
        if (nsd == 3)
        {
          utl::getAttribute(child,"direction",load.ldof);
          IFEM::cout <<" direction = "<< (int)load.ldof;
        }
        IFEM::cout <<" load = "<< load.pload << std::endl;
        myLoads.push_back(load);
      }
    }

    else if (!strcasecmp(child->Value(),"pressure")) {
      std::string set, type;
      utl::getAttribute(child,"set",set);
      int code = this->getUniquePropertyCode(set,1);
      if (code == 0) utl::getAttribute(child,"code",code);

      if (child->FirstChild() && code > 0) {
        utl::getAttribute(child,"type",type,true);
        IFEM::cout <<"\tPressure code "<< code;
        if (!type.empty()) IFEM::cout <<" ("<< type <<")";
        myScalars[code] = utl::parseRealFunc(child->FirstChild()->Value(),type);
        this->setPropertyType(code,Property::BODYLOAD);
        IFEM::cout << std::endl;
      }
    }

    else if (!strcasecmp(child->Value(),"anasol")) {
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


bool SIMLinElKL::initMaterial (size_t propInd)
{
  if (propInd >= mVec.size()) propInd = mVec.size()-1;

  KirchhoffLove* klp = dynamic_cast<KirchhoffLove*>(myProblem);
  if (!klp) return false;

  klp->setMaterial(mVec[propInd]);
  if (tVec[propInd] != 0.0)
    klp->setThickness(tVec[propInd]);

  return true;
}


bool SIMLinElKL::initBodyLoad (size_t patchInd)
{
  KirchhoffLove* klp = dynamic_cast<KirchhoffLove*>(myProblem);
  if (!klp) return false;

  SclFuncMap::const_iterator it = myScalars.find(0);
  for (const Property& prop : myProps)
    if (prop.pcode == Property::BODYLOAD && prop.patch == patchInd)
      if ((it = myScalars.find(prop.pindx)) != myScalars.end()) break;

  klp->setPressure(it == myScalars.end() ? 0 : it->second);
  return true;
}


bool SIMLinElKL::initNeumann (size_t propInd)
{
  KirchhoffLove* klp = dynamic_cast<KirchhoffLove*>(myProblem);
  if (!klp) return false;

  return true;
}


void SIMLinElKL::preprocessA ()
{
  if (!myProblem)
  {
    if (nsd == 3)
      myProblem = new KirchhoffLoveShell();
    else
      myProblem = new KirchhoffLovePlate();
  }

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


bool SIMLinElKL::preprocessB ()
{
  // Preprocess the nodal point loads
  for (PloadVec::iterator p = myLoads.begin(); p != myLoads.end();)
    if ((p->inod = this->evalPoint(p->xi,p->X,nullptr,p->patch,true)) < 1)
    {
      p = myLoads.erase(p);
      std::cerr <<"  ** SIMLinElKL::preprocess: Load point ("
                << p->xi[0] <<','<< p->xi[1] <<") on patch #"<< p->patch
                <<" is not a nodal point (ignored)."<< std::endl;
    }
    else
    {
      int ipt = 1 + (int)(p-myLoads.begin());
      if (ipt == 1) IFEM::cout <<'\n';
      IFEM::cout <<"Load point #"<< ipt <<": patch #"<< p->patch
                 <<" (u,v)=("<< p->xi[0] <<','<< p->xi[1]
                 <<"), node #"<< p->inod <<", X = "<< p->X;
      if (nsd == 3)
        IFEM::cout <<", direction = "<< (int)p->ldof;
      IFEM::cout << std::endl;
      ++p;
    }

  return true;
}


bool SIMLinElKL::assembleDiscreteTerms (const IntegrandBase*, const TimeDomain&)
{
  SystemVector* b = myEqSys->getVector();
  if (!b) return false;

  for (const PointLoad& load : myLoads)
  {
    double pload[3] = { 0.0, 0.0, 0.0 };
    pload[load.ldof-1] = load.pload;
    if (!mySam->assembleSystem(*b,pload,load.inod))
      return false;
  }

  return true;
}


double SIMLinElKL::externalEnergy (const Vectors& psol) const
{
  double energy = this->SIMbase::externalEnergy(psol);

  // External energy from the nodal point loads
  const int* madof = mySam->getMADOF();
  for (const PointLoad& load : myLoads)
    energy += load.pload * psol.front()(madof[load.inod-1] + load.ldof-1);

  return energy;
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
