// $Id$
//==============================================================================
//!
//! \file SIMKLShell.C
//!
//! \date Sep 16 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based FE analysis of Kirchhoff-Love shells.
//!
//==============================================================================

#include "SIMKLShell.h"
#include "KirchhoffLoveShell.h"
#include "LinIsotropic.h"
#include "AlgEqSystem.h"
#include "ASMbase.h"
#include "SAM.h"
#include "Functions.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "Property.h"
#include "tinyxml.h"


SIMKLShell::SIMKLShell (bool shell)
{
  if (shell) nsd = 3;

  nf[0] = shell ? 3 : 1;
  aCode[0] = aCode[1] = aCode[2] = 0;
}


SIMKLShell::~SIMKLShell ()
{
  for (PointLoad& load : myLoads)
    delete load.p;

  // To prevent the SIMbase destructor try to delete already deleted functions
  for (int i = 0; i < 3; i++)
    if (aCode[i] > 0) myScalars.erase(aCode[i]);
}


KirchhoffLove* SIMKLShell::getProblem (int)
{
  KirchhoffLove* klp = dynamic_cast<KirchhoffLove*>(myProblem);
  if (!klp)
    myProblem = klp = new KirchhoffLoveShell();

  return klp;
}


Elasticity* SIMKLShell::getIntegrand ()
{
  return dynamic_cast<Elasticity*>(myProblem);
}


void SIMKLShell::clearProperties()
{
  for (PointLoad& load : myLoads)
    delete load.p;

  myLoads.clear();

  this->SIMElasticity<SIM2D>::clearProperties();
}


bool SIMKLShell::parse (char* keyWord, std::istream& is)
{
  KirchhoffLove* klp = this->getProblem();

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
      myLoads[i].p = new ConstantFunc(atof(strtok(NULL," ")));
      IFEM::cout <<"\n\tPoint "<< i+1 <<": P"<< myLoads[i].patch
                 <<" xi = "<< myLoads[i].xi[0] <<" "<< myLoads[i].xi[1];
      if (nsd == 3)
      {
        if ((cline = strtok(NULL," ")))
          myLoads[i].ldof.second = atoi(cline);
        IFEM::cout <<" direction = "<< myLoads[i].ldof.second;
      }
      IFEM::cout <<" load = "<< (*myLoads[i].p)(0.0);
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

  else
    return this->SIM2D::parse(keyWord,is);

  return true;
}


bool SIMKLShell::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"kirchhofflove"))
    return this->SIM2D::parse(elem);

  int version = 1;
  utl::getAttribute(elem,"version",version);
  KirchhoffLove* klp = this->getProblem(version);

  bool ok = true;
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child && ok; child = child->NextSiblingElement())

    if (!strcasecmp(child->Value(),"gravity"))
    {
      double g = 0.0;
      utl::getAttribute(child,"g",g);
      klp->setGravity(g);
      IFEM::cout <<"\tGravitation constant: "<< g << std::endl;
    }

    else if (!strcasecmp(child->Value(),"isotropic"))
    {
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

    else if (!strcasecmp(child->Value(),"pointload") && child->FirstChild())
    {
      myLoads.resize(myLoads.size()+1);
      PointLoad& load = myLoads.back();
      utl::getAttribute(child,"patch",load.patch);
      utl::getAttribute(child,"xi",load.xi[0]);
      utl::getAttribute(child,"eta",load.xi[1]);
      IFEM::cout <<"\tPoint: P"<< load.patch
                 <<" xi = "<< load.xi[0] <<" "<< load.xi[1];
      if (nsd == 3)
      {
        utl::getAttribute(child,"direction",load.ldof.second);
        IFEM::cout <<" direction = "<< load.ldof.second;
      }

      std::string type("constant");
      utl::getAttribute(child,"type",type);
      if (type == "constant")
      {
        load.p = new ConstantFunc(atof(child->FirstChild()->Value()));
        IFEM::cout <<" load = "<< (*load.p)(0.0) << std::endl;
      }
      else
      {
        IFEM::cout <<" Load: ";
        load.p = utl::parseTimeFunc(child->FirstChild()->Value(),type);
      }

      // Negate the local DOF flag to signal that element loads are allowed
      bool allowElementPointLoad = false;
      utl::getAttribute(child,"onElement",allowElementPointLoad);
      if (allowElementPointLoad) load.ldof.second *= -1;
    }

    else if (!strcasecmp(child->Value(),"pressure") && child->FirstChild())
    {
      std::string set, type;
      utl::getAttribute(child,"set",set);
      int code = this->getUniquePropertyCode(set,1);
      if (code == 0) utl::getAttribute(child,"code",code);
      if (code > 0)
      {
        utl::getAttribute(child,"type",type,true);
        IFEM::cout <<"\tPressure code "<< code;
        if (!type.empty()) IFEM::cout <<" ("<< type <<")";
        myScalars[code] = utl::parseRealFunc(child->FirstChild()->Value(),type);
        this->setPropertyType(code,Property::BODYLOAD);
        IFEM::cout << std::endl;
        if (!klp->haveLoads('I'))
          klp->setPressure(myScalars[code]);
      }
    }

    else
      ok = this->SIM2D::parse(child);

  return ok;
}


bool SIMKLShell::initMaterial (size_t propInd)
{
  if (propInd >= mVec.size()) propInd = mVec.size()-1;

  KirchhoffLove* klp = dynamic_cast<KirchhoffLove*>(myProblem);
  if (!klp) return false;

  klp->setMaterial(mVec[propInd]);
  if (tVec[propInd] != 0.0)
    klp->setThickness(tVec[propInd]);

  return true;
}


bool SIMKLShell::initBodyLoad (size_t patchInd)
{
  KirchhoffLove* klp = dynamic_cast<KirchhoffLove*>(myProblem);
  if (!klp) return false;

  klp->setPressure();
  SclFuncMap::const_iterator it = myScalars.find(0);
  if (it != myScalars.end()) klp->setPressure(it->second);

  for (const Property& prop : myProps)
    if (prop.pcode == Property::BODYLOAD && prop.patch == patchInd)
      if ((it = myScalars.find(prop.pindx)) != myScalars.end())
        if (it->second) klp->setPressure(it->second);

  return true;
}


bool SIMKLShell::initNeumann (size_t propInd)
{
  KirchhoffLove* klp = dynamic_cast<KirchhoffLove*>(myProblem);
  if (!klp) return false;

  VecFuncMap::const_iterator vit = myVectors.find(propInd);
  TracFuncMap::const_iterator tit = myTracs.find(propInd);

  if (vit != myVectors.end())
    klp->setTraction(vit->second);
  else if (tit != myTracs.end())
    klp->setTraction(tit->second);
  else
    return false;

  return true;
}


void SIMKLShell::preprocessA ()
{
  this->getProblem();
  this->printProblem();
}


bool SIMKLShell::preprocessB ()
{
  // Preprocess the nodal point loads, if any
  if (myLoads.empty())
    return true;

  IFEM::cout <<'\n';
  bool ok = true;
  int ipt = 0;
  for (PointLoad& pl : myLoads)
  {
    double prms[2];
    int iclose = 0;
    int imatch = this->evalPoint(pl.xi,pl.X,prms,pl.patch,true);
    if (imatch >= 0 && pl.ldof.second < 0)
      pl.ldof.first = this->findElementContaining(prms,pl.patch);
    else if (imatch == 0 && (iclose = this->findClosestNode(pl.X)) > 0)
      pl.ldof.first = iclose;
    else
      pl.ldof.first = imatch;

    if (pl.ldof.first <= 0)
    {
      std::cerr <<" *** SIMKLShell::preprocessB: Load point ("
                << pl.xi[0] <<','<< pl.xi[1] <<") on patch #"<< pl.patch
                <<" is not a nodal point."<< std::endl;
      ok = false;
      continue;
    }

    memcpy(pl.xi,prms,2*sizeof(double));
    IFEM::cout <<"Load point #"<< ++ipt <<": patch #"<< pl.patch
               <<" (u,v)=("<< pl.xi[0] <<','<< pl.xi[1];
    if (pl.ldof.second < 0)
      IFEM::cout <<") on element #";
    else if (iclose > 0)
      IFEM::cout <<"), (closest) node #";
    else
      IFEM::cout <<"), node #";
    IFEM::cout << pl.ldof.first <<", X = "<< pl.X;
    if (iclose > 0)
      IFEM::cout <<" (Xnod = "<< this->getNodeCoord(iclose) <<")";
    if (nsd == 3)
      IFEM::cout <<", direction = "<< abs(pl.ldof.second);
    IFEM::cout << std::endl;
  }

  return ok;
}


bool SIMKLShell::assembleDiscreteTerms (const IntegrandBase*,
                                        const TimeDomain& time)
{
  SystemVector* b = myEqSys->getVector();
  if (!b) return false;

  bool ok = true;
  for (const PointLoad& load : myLoads)
    if (load.ldof.second > 0)
      ok &= mySam->assembleSystem(*b,(*load.p)(time.t),load.ldof);
    else // This is an element point load
      ok &= this->assemblePoint(load.patch,load.xi,(*load.p)(time.t),
                                -load.ldof.second);

  if (ok && time.first && time.it == 0)
  {
    Vector extLoad;
    if (mySam->expandSolution(*b,extLoad,0.0))
    {
      std::streamsize oldPrec = IFEM::cout.precision(15);
      IFEM::cout <<"  * Sum external load:";
      for (unsigned char d = 0; d < nf[0]; d++)
        IFEM::cout <<" "<< extLoad.sum(d,nf[0]);
      IFEM::cout << std::endl;
      IFEM::cout.precision(oldPrec);
    }
  }

  return ok;
}


double SIMKLShell::externalEnergy (const Vectors& u,
                                   const TimeDomain& time) const
{
  double energy = this->SIMbase::externalEnergy(u,time);

  // External energy from the nodal point loads
  const int* madof = mySam->getMADOF();
  for (const PointLoad& load : myLoads)
    if (load.ldof.second > 0)
    {
      int idof = madof[load.ldof.first-1] + load.ldof.second-1;
      energy += (*load.p)(time.t) * u.front()(idof);
    }
    else if (load.ldof.second < 0) // This is an element point load
    {
      Vector v = this->SIMgeneric::getSolution(u.front(),load.xi,0,load.patch);
      if (-load.ldof.second <= (int)v.size())
        energy += (*load.p)(time.t) * v(-load.ldof.second);
    }

  return energy;
}


bool SIMKLShell::assemblePoint (int patch, const double* u, double f, int ldof)
{
  ASMbase* pch = this->getPatch(patch,true);
  if (!pch) return false;

  Vec3 Fvec; Fvec(ldof) = f;
  this->setMode(SIM::RHS_ONLY);
  return pch->diracPoint(*myProblem,*myEqSys,u,Fvec);
}
