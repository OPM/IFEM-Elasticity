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
#include "Elasticity.h"
#include "ElasticityUtils.h"
#include "KirchhoffLoveShell.h"
#include "LinIsotropic.h"
#include "AlgEqSystem.h"
#include "ASMs2D.h"
#ifdef HAS_LRSPLINE
#include "ASMu2D.h"
#endif
#include "SAM.h"
#include "Functions.h"
#include "IFEM.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "Property.h"
#include "tinyxml.h"


SIMKLShell::SIMKLShell (const char* heading, bool isShell)
{
  if (heading) myHeading = heading;
  if (isShell) nsd = 3;

  myContext = "kirchhofflove";

  nf[0] = isShell ? 3 : 1;
  if (opt.discretization == ASM::LRSpline)
  {
    nf.push_back('C');
    nf.push_back('1');
  }

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
  // To prevent SIMbase::clearProperties deleting the analytical solution
  for (int i = 0; i < 3; i++)
    if (aCode[i] > 0) myScalars.erase(aCode[i]);

  aCode[0] = aCode[1] = aCode[2] = 0;

  KirchhoffLove* klp = dynamic_cast<KirchhoffLove*>(myProblem);
  if (klp)
  {
    klp->setMaterial(nullptr);
    klp->setPressure(nullptr);
  }

  for (Material* mat : mVec)
    delete mat;
  for (PointLoad& load : myLoads)
    delete load.p;

  tVec.clear();
  mVec.clear();
  myLoads.clear();

  this->SIM2D::clearProperties();
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

  else if (!strncasecmp(keyWord,"ANASOL",6))
    return this->parseAnaSol(keyWord,is);

  else
    return this->SIM2D::parse(keyWord,is);

  return true;
}


bool SIMKLShell::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),myContext.c_str()))
    return this->SIMElasticity<SIM2D>::parse(elem);

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

    else if (!strcasecmp(child->Value(),"lineload") && child->FirstChild())
    {
      std::string set, type;
      utl::getAttribute(child,"set",set);
      int code = this->getUniquePropertyCode(set,1);
      if (code == 0) utl::getAttribute(child,"code",code);
      if (code > 0)
      {
        utl::getAttribute(child,"type",type,true);
        IFEM::cout <<"\tLine load code "<< code;
        if (!type.empty()) IFEM::cout <<" ("<< type <<")";
        myScalars[code] = utl::parseRealFunc(child->FirstChild()->Value(),type);
        this->setPropertyType(code,Property::OTHER);
        klp->setLineLoad(myScalars[code]);
        if (utl::getAttribute(child,"u",lineLoad.u))
        {
          lineLoad.direction = 2;
          utl::getAttribute(child,"v0",lineLoad.range.first);
          utl::getAttribute(child,"v1",lineLoad.range.second);
          IFEM::cout <<"\n\tat u="<< lineLoad.u
                     <<", v in ["<< lineLoad.range.first
                     <<","<< lineLoad.range.second <<"]";
        }
        else if (utl::getAttribute(child,"v",lineLoad.u))
        {
          lineLoad.direction = 1;
          utl::getAttribute(child,"u0",lineLoad.range.first);
          utl::getAttribute(child,"u1",lineLoad.range.second);
          IFEM::cout <<"\n\tat v="<< lineLoad.u
                     <<", u in ["<< lineLoad.range.first
                     <<","<< lineLoad.range.second <<"]";
        }
        IFEM::cout << std::endl;
      }
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

    else if (!strcasecmp(child->Value(),"rigid"))
      ok &= this->parseRigid(child,this);

    else if (!strcasecmp(child->Value(),"anasol"))
      ok = this->parseAnaSol(child);

    else if (!klp->parse(child))
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
  // Check if we have neumann conditions (for Gauss point visualization).
  // Need to do this to enable calculation of the number of boundary
  // integration points before any traction field is assigned.
  VecFunc* dummy = nullptr;
  for (const Property& prop : myProps)
    if (prop.pcode == Property::NEUMANN)
      ++dummy;

  this->getProblem()->setTraction(dummy);
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


bool SIMKLShell::assembleDiscreteTerms (const IntegrandBase* itg,
                                        const TimeDomain& time)
{
  if (!myEqSys) return true; // silently ignore if no equation system defined

  SystemVector* b = myEqSys->getVector();
  if (!b) return false;

  bool ok = true;
  if (itg != myProblem)
    return ok;

  SIM::SolutionMode mode = itg->getMode();
  if (!myLoads.empty())
    this->setMode(SIM::RHS_ONLY);

  for (const PointLoad& load : myLoads)
    if (load.ldof.second > 0)
      ok &= mySam->assembleSystem(*b,(*load.p)(time.t),load.ldof);
    else // This is an element point load
      ok &= this->assemblePoint(load.patch,load.xi,(*load.p)(time.t),
                                -load.ldof.second);

  if (ok && time.first && time.it == 0 && mode != SIM::ARCLEN)
  {
    Vector extLoad;
    if (mySam->expandSolution(*b,extLoad,0.0))
    {
      std::streamsize oldPrec = IFEM::cout.precision(15);
      IFEM::cout <<"   * Sum external load:";
      for (unsigned char d = 0; d < nf[0]; d++)
        IFEM::cout <<" "<< extLoad.sum(d,nf[0]);
      IFEM::cout << std::endl;
      IFEM::cout.precision(oldPrec);
    }
  }

  if (!myLoads.empty() && mode == SIM::ARCLEN)
    b = myEqSys->getVector(1); // External load gradient for arc-length driver
  else
    b = nullptr;

  if (b)
  {
    // Assemble external nodal point load gradient at current time step
    this->setMode(mode);
    static_cast<KirchhoffLove*>(myProblem)->setLoadGradientMode();
    for (const PointLoad& load : myLoads)
      if (load.ldof.second > 0)
        ok &= mySam->assembleSystem(*b,load.p->deriv(time.t),load.ldof);
      else // This is an element point load
        ok &= this->assemblePoint(load.patch,load.xi,load.p->deriv(time.t),
                                  -load.ldof.second);
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


void SIMKLShell::shiftGlobalNums (int nshift, int)
{
  for (PointLoad& load : myLoads)
    if (load.ldof.second > 0)
      load.ldof.first += nshift;
}


bool SIMKLShell::assemblePoint (int patch, const double* u, double f, int ldof)
{
  ASMbase* pch = this->getPatch(patch,true);
  if (!pch) return false;

  Vec3 Fvec; Fvec(ldof) = f;
  return pch->diracPoint(*myProblem,*myEqSys,u,Fvec);
}


bool SIMKLShell::getExtLoad (RealArray& extloa, const TimeDomain& time) const
{
  extloa.resize(nf[0]);
  for (size_t i = 0; i < nf[0]; i++)
    extloa[i] = this->extractScalar(i);

  for (const PointLoad& load : myLoads)
    if (load.ldof.second > 0 && load.ldof.second < nf[0])
      extloa[load.ldof.second-1] += (*load.p)(time.t);

  return true;
}


void SIMKLShell::printStep (int istep, const TimeDomain& time) const
{
  adm.cout <<"\n  step="<< istep <<"  time="<< time.t;

  RealArray extLo;
  if (myProblem->getMode() == SIM::ARCLEN && this->getExtLoad(extLo,time))
  {
    adm.cout <<"  Sum(Fex) =";
    for (size_t d = 0; d < extLo.size(); d++)
      adm.cout <<" "<< utl::trunc(extLo[d]);
  }

  adm.cout << std::endl;
}


void SIMKLShell::printNormGroup (const Vector& norm, const Vector& rNorm,
                                 const std::string& prjName) const

{
  if (!dynamic_cast<KirchhoffLoveShell*>(myProblem))
  {
    Elastic::printNorms(norm,rNorm,prjName,this);
    return;
  }

  if (utl::trunc(norm(1)) == 0.0)
    return;

  // Special print for shell problems (no analytical solution)
  IFEM::cout <<"\n\n>>> Error estimates based on "<< prjName <<" <<<"
             <<"\nEnergy norm |u^r| = a(u^r,u^r)^0.5   : "<< norm(1)
             <<"\nError norm a(e,e)^0.5, e=u^r-u^h     : "<< norm(2)
             <<"\n- relative error (% of |u^r|) : "<< 100.0*norm(2)/norm(1);
  if (utl::trunc(norm(3)) > 0.0)
    IFEM::cout <<"\nL2 norm |n^r| = (n^r,n^r)^0.5        : "<< norm(3)
               <<"\nL2 error (e,e)^0.5, e=n^r-n^h        : "<< norm(5)
               <<"\n- relative error (% of |n^r|) : "<< 100.0*norm(5)/norm(3);
  if (utl::trunc(norm(4)) > 0.0)
    IFEM::cout <<"\nL2 norm |m^r| = (m^r,m^r)^0.5        : "<< norm(4)
               <<"\nL2 error (e,e)^0.5, e=m^r-m^h        : "<< norm(6)
               <<"\n- relative error (% of |m^r|) : "<< 100.0*norm(6)/norm(4);
}


ASM::InterfaceChecker* SIMKLShell::getInterfaceChecker (size_t pidx) const
{
  ASM::InterfaceChecker* ichk = nullptr;

  const ASMs2D* spch = dynamic_cast<const ASMs2D*>(this->getPatch(pidx+1));
  if (spch)
    ichk = new LineLoadChecker<ASMs2D>(*spch);
#ifdef HAS_LRSPLINE
  else
  {
    const ASMu2D* upch = dynamic_cast<const ASMu2D*>(this->getPatch(pidx+1));
    if (upch)
      ichk = new LineLoadChecker<ASMu2D>(*upch);
  }
#endif

  if (ichk)
  {
    if (lineLoad.direction == 1)
      dynamic_cast<ElmBorderChk*>(ichk)->setDomainEW(lineLoad.u,lineLoad.range);
    else if (lineLoad.direction == 2)
      dynamic_cast<ElmBorderChk*>(ichk)->setDomainNS(lineLoad.u,lineLoad.range);
  }

  return ichk;
}


ElmBorderChk::ElmBorderChk ()
{
  umin = umax = 0.5;
  vmin = 0.25;
  vmax = 0.75;
}


void ElmBorderChk::setDomainNS (double u, const Doubles& rng)
{
  umin = umax = u;
  vmin = rng.first;
  vmax = rng.second;
}


void ElmBorderChk::setDomainEW (double v, const Doubles& rng)
{
  vmin = vmax = v;
  umin = rng.first;
  umax = rng.second;
}


short int ElmBorderChk::maskBorder (double u0, double u1,
                                    double v0, double v1) const
{
  if (umin == umax && v0 >= vmin && v1 <= vmax)
  {
    if (u0 == umin)
      return 1;
    else if (u1 == umax)
      return 2;
  }
  else if (vmin == vmax && u0 >= umin && u1 <= umax)
  {
    if (v0 == vmin)
      return 4;
    else if (v1 == vmax)
      return 8;
  }

  return 0;
}
