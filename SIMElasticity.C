// $Id$
//==============================================================================
//!
//! \file SIMElasticity.C
//!
//! \date Dec 04 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#include "SIMElasticity.h"

#include "Elasticity.h"
#include "ElasticityUtils.h"
#include "SIMgeneric.h"
#include "MaterialBase.h"

#include "AnaSol.h"
#include "ASMbase.h"
#include "ForceIntegrator.h"
#include "Functions.h"
#include "IFEM.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "TractionField.h"
#include "TimeStep.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "VTF.h"

#include "tinyxml2.h"


//! Plane strain/stress option for 2D problems.
bool Elastic::planeStrain = false;
//! Axisymmtry option for 2D problems.
bool Elastic::axiSymmetry = false;
//! Option for Gauss point output to VTF.
bool Elastic::GIpointsVTF = false;
//! Time for function evaluation in linear problems.
double Elastic::time = 1.0;


void Elastic::printNorms (const Vector& gNorm, const Vector& rNorm,
                          const std::string& name, const SIMbase* model)
{
  const char* uRef = model->haveAnaSol() ? "|u|)   : " : "|u^r|) : ";
  double Rel = 100.0/(model->haveAnaSol() ? rNorm(3) : gNorm(1));

  IFEM::cout <<"\n\n>>> Error estimates based on "<< name <<" <<<"
             <<"\nEnergy norm |u^r| = a(u^r,u^r)^0.5   : "<< gNorm(1)
             <<"\nError norm a(e,e)^0.5, e=u^r-u^h     : "<< gNorm(2)
             <<"\n- relative error (% of "<< uRef << gNorm(2)*Rel;

  bool haveResErr = gNorm.size() >= (model->haveAnaSol() ? 8 : 5);
  if (haveResErr)
  {
    IFEM::cout <<"\nResidual error (r(u^r) + J(u^r))^0.5 : "<< gNorm(5)
               <<"\n- relative error (% of "<< uRef << gNorm(5)*Rel;
    if (gNorm.size() >= 9)
      IFEM::cout <<"\nJump term J(u^r)^0.5          : "<< gNorm(6)
                 <<"\n- relative error (% of "<< uRef << gNorm(6)*Rel;
  }

  if (model->haveAnaSol())
  {
    double exaErr = gNorm(gNorm.size() - (haveResErr ? 2 : 1));
    IFEM::cout <<"\nExact error a(e,e)^0.5, e=u-u^r      : "<< exaErr
               <<"\n- relative error (% of "<< uRef << exaErr*Rel
               <<"\nEffectivity index             : "<< gNorm(2)/rNorm(4);
    if (haveResErr)
      IFEM::cout <<"\nEffectivity index, theta^EX          : "
                 << (gNorm(2)+exaErr)/rNorm(4)
                 <<"\nEffectivity index, theta^RES         : "
                 << (gNorm(2)+gNorm(5))/rNorm(4);
  }
}


void Elastic::printBoundaryForces (const Vector& sf, RealArray& weights,
                                   const std::map<int,Vec3>& bCode,
                                   const SIMgeneric* model,
                                   bool indented)
{
  if (bCode.size() > 1 && weights.empty())
  {
    std::vector<int> glbNodes;
    weights.resize(model->getNoNodes());
    for (const std::pair<const int,Vec3>& c : bCode)
    {
      model->getBoundaryNodes(c.first,glbNodes);
      for (int inod : glbNodes) ++weights[inod-1];
    }
    for (double& w : weights)
      if (w > 1.0) w = 1.0/w;
  }

  int isec = 0;
  for (const std::pair<const int,Vec3>& c : bCode)
  {
    ++isec;
    RealArray force = model->getInterfaceForces(sf,weights,c.first);
    if (Vector(force).normInf() > utl::zero_print_tol)
    {
      IFEM::cout << (indented ? "\n  " : "\n")
                 <<"Interface force at section "<< isec <<":";
      for (double f : force) IFEM::cout <<" "<< utl::trunc(f);
    }
  }

  IFEM::cout << std::endl;
}


template<class Dim>
SIMElasticity<Dim>::SIMElasticity (bool checkRHS) : Dim(Dim::dimension,checkRHS)
{
  myContext = "elasticity";
  plotRgd = printed = false;
  aCode = 0;
}


template<class Dim>
SIMElasticity<Dim>::~SIMElasticity ()
{
  // To prevent the SIMbase destructor try to delete already deleted functions
  if (aCode > 0)
    Dim::myVectors.erase(aCode);

  for (Material* mat : mVec)
    delete mat;
}


template<class Dim>
bool SIMElasticity<Dim>::printProblem () const
{
  // Avoid printing problem definition more than once
  if (printed) return false;

  return printed = this->Dim::printProblem();
}


/*!
  This method is reimplemented to print out the external load in the beginning
  of the load step in case of arc-length solution driver.
*/

template<class Dim>
void SIMElasticity<Dim>::printStep (int istep, const TimeDomain& time) const
{
  Dim::adm.cout <<"\n  step="<< istep <<"  time="<< time.t;

  if (Dim::myProblem->getMode() == SIM::ARCLEN)
  {
    RealArray extLo;
    if (this->extractScalars(extLo))
    {
      Dim::adm.cout <<"  Sum(Fex) =";
      for (double f : extLo) Dim::adm.cout <<" "<< utl::trunc(f);
    }
  }

  Dim::adm.cout << std::endl;
}


template<class Dim>
bool SIMElasticity<Dim>::advanceStep (TimeStep& tp)
{
  Elasticity* elp = dynamic_cast<Elasticity*>(Dim::myProblem);
  if (elp)
    elp->advanceStep(tp.time.dt,tp.time.dtn);

  return true;
}


template<class Dim>
void SIMElasticity<Dim>::clearProperties ()
{
  // To prevent SIMbase::clearProperties deleting the analytical solution
  if (aCode > 0)
    Dim::myVectors.erase(aCode);
  aCode = 0;

  bCode.clear();

  Elasticity* elp = dynamic_cast<Elasticity*>(Dim::myProblem);
  if (elp)
  {
    elp->setMaterial(nullptr);
    elp->setBodyForce(nullptr);
    elp->setTraction((VecFunc*)nullptr);
    elp->setTraction((TractionFunc*)nullptr);
    elp->addExtrFunction(nullptr);
  }

  for (Material* mat : mVec)
    delete mat;
  mVec.clear();

  this->Dim::clearProperties();
}


/*!
  The boundaries for which the surface traction resultants are calculated
  are identified by the property set codes in \ref bCode, which are assigned
  values by parsing the `<boundaryforce>` tags in the input file.
*/

template<class Dim>
bool SIMElasticity<Dim>::calcBouForces (Real2DMat& f, const Vectors& sol)
{
  f.clear();
  f.reserve(bCode.size());
  TimeDomain time;
  for (const std::pair<const int,Vec3>& c : bCode)
  {
    f.push_back(SIM::getBoundaryForce(sol,this,c.first,time,&c.second));
    Dim::adm.allReduceAsSum(f.back());
  }

  return !f.empty();
}


/*!
  The boundary for which the surface traction resultant is calculated
  is identified by the first property set code in \ref bCode, which is assigned
  value by parsing the first `<boundaryforce>` tag in the input file.
*/

template<class Dim>
bool SIMElasticity<Dim>::getBoundaryForce (RealArray& f,
                                           const Vectors& sol,
                                           const TimeStep& tp)
{
  if (bCode.empty())
    return false;

  f = SIM::getBoundaryForce(sol,this,bCode.begin()->first,tp.time);
  Dim::adm.allReduceAsSum(f);
  return true;
}


/*!
  The boundaries for which the reaction forces are returned
  are identified by the property set codes in \ref bCode, which are assigned
  values by parsing the `<boundaryforce>` tags in the input file.
*/

template<class Dim>
bool SIMElasticity<Dim>::getBoundaryReactions (Real2DMat& rf)
{
  rf.resize(bCode.size());
  if (bCode.empty())
    return false;

  size_t i = 0;
  bool ok = true;
  for (const std::pair<const int,Vec3>& c : bCode)
  {
    ok &= this->getCurrentReactions(rf[i],{},c.first);
    Dim::adm.allReduceAsSum(rf[i++]);
  }

  return ok;
}


template<class Dim>
bool SIMElasticity<Dim>::getBoundaryReactions (Vector& rf, size_t bindex)
{
  Real2DMat rtmp;
  if (!this->getBoundaryReactions(rtmp) || bindex > rtmp.size())
    return false;
  else if (bindex > 0)
    rf = rtmp[bindex];
  else
  {
    rf = rtmp.front();
    for (size_t i = 1; i < rtmp.size(); i++)
      rf.add(rtmp[i]);
  }

  return true;
}


template<class Dim>
bool SIMElasticity<Dim>::haveBoundaryReactions (bool reactionsOnly) const
{
  for (const std::pair<const int,Vec3>& c : bCode)
    if (!reactionsOnly || this->haveReactions(c.first))
      return true;

  return false;
}


template<class Dim>
bool SIMElasticity<Dim>::haveAnaSol () const
{
  return (Dim::mySol && Dim::mySol->getStressSol());
}


/*!
  This method is reimplemented inserting a call to the method getIntegrand().
  This makes sure the integrand has been allocated in case of minimum input.
  It also resolves inhomogeneous boundary condition fields in case they are
  derived from the analytical solution.
*/

template<class Dim>
void SIMElasticity<Dim>::preprocessA ()
{
  ElasticBase* elInt = this->getIntegrand();

  this->printProblem();

  // Deactivate principal stress ouput for Lagrange/Spectral interpolations
  if (Dim::opt.discretization < ASM::Spline)
    Elasticity::wantPrincipalStress = false;

  if (Dim::dualField)
    static_cast<Elasticity*>(elInt)->setDualRHS(Dim::dualField);

  if (!Dim::mySol) return;

  // Define analytical boundary condition fields
  for (Property& p : Dim::myProps)
    if (p.pcode == Property::DIRICHLET_ANASOL)
    {
      VecFunc* vecField = Dim::mySol->getVectorSol();
      if (!vecField)
        p.pcode = Property::UNDEFINED;
      else if (aCode == abs(p.pindx))
        p.pcode = Property::DIRICHLET_INHOM;
      else if (aCode == 0)
      {
        aCode = abs(p.pindx);
        Dim::myVectors[aCode] = vecField;
        p.pcode = Property::DIRICHLET_INHOM;
      }
      else
        p.pcode = Property::UNDEFINED;
    }
    else if (p.pcode == Property::NEUMANN_ANASOL)
    {
      STensorFunc* stressField = Dim::mySol->getStressSol();
      if (stressField)
      {
        p.pcode = Property::NEUMANN;
        if (Dim::myTracs.find(p.pindx) == Dim::myTracs.end())
          Dim::myTracs[p.pindx] = new TractionField(*stressField);
      }
      else
        p.pcode = Property::UNDEFINED;
    }
}


/*!
  This method creates the multi-point constraint equations representing the
  rigid couplings in the model, if any.
*/

template<class Dim>
bool SIMElasticity<Dim>::preprocessBeforeAsmInit (int& ngnod)
{
  return this->addRigidMPCs(this,ngnod);
}


/*!
  This method is reimplemented to ensure that threading groups are established
  for the patch faces subjected to boundary force integration.
  In addition, the reference point for moment calculation \b X0 of each boundary
  is calculated based on the control/nodal point coordinates.
*/

template<class Dim>
bool SIMElasticity<Dim>::preprocessB ()
{
  size_t iSec = 0;
  for (std::pair<const int,Vec3>& code : bCode)
  {
    for (const Property& p : Dim::myProps)
      if (p.pindx == code.first)
        this->generateThreadGroups(p,Dim::msgLevel < 2);

    IntVec nodes;
    Vec3Vec Xnodes;
    this->getBoundaryNodes(code.first,nodes,&Xnodes);
    if (nodes.empty()) continue;

    // Find the centre of all boundary control/nodal points
    Vec3& X0 = code.second;
    for (const Vec3& X : Xnodes) X0 += X;
    X0 /= Xnodes.size();

    // Find the location of the point which is furthest away from the centre
    Vec3 X1(X0);
    double d, dmax = 0.0;
    for (const Vec3& X : Xnodes)
      if ((d = (X-X0).length2()) > dmax)
      {
        dmax = d;
        X1 = X;
      }

    // Find the location of the point which is furthest away from X1
    Vec3 X2(X1);
    for (const Vec3& X : Xnodes)
      if ((d = (X-X1).length2()) > dmax)
      {
        dmax = d;
        X2 = X;
      }

    // Assuming X1 and X2 now are the end points of the straight 1D boundary,
    // or the diameter of the smallest subscribing circle of a 2D boundary,
    // the reference point X0 is taken as the mid-point between X1 and X2
    X0 = 0.5*(X1+X2);
    IFEM::cout <<"Boundary section "<< ++iSec <<": X0 = "<< X0 << std::endl;
  }

  return true;
}


template<class Dim>
bool SIMElasticity<Dim>::parseAnaSol (char*, std::istream&)
{
  std::cerr <<" *** SIMElasticity::parse: No analytical solution available."
            << std::endl;
  return false;
}


template<class Dim>
bool SIMElasticity<Dim>::parseAnaSol (const tinyxml2::XMLElement*)
{
  std::cerr <<" *** SIMElasticity::parse: No analytical solution available."
            << std::endl;
  return false;
}


template<class Dim>
bool SIMElasticity<Dim>::parse (char* keyWord, std::istream& is)
{
  char* cline = nullptr;
  int nmat = 0;
  int nConstPress = 0;
  int nLinearPress = 0;

  if (!strncasecmp(keyWord,"ISOTROPIC",9))
  {
    nmat = atoi(keyWord+10);
    IFEM::cout <<"\nNumber of isotropic materials: "<< nmat << std::endl;
    ElasticBase* elInt = this->getIntegrand();
    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      IFEM::cout <<"\tMaterial code "<< code <<": ";
      if (code > 0)
        this->setPropertyType(code,Property::MATERIAL,mVec.size());
      bool planeStrain = Dim::dimension == 2 ? Elastic::planeStrain : true;
      mVec.push_back(elInt->parseMatProp((char*)nullptr,planeStrain));
      IFEM::cout << std::endl;
    }
  }

  else if (!strncasecmp(keyWord,"GRAVITY",7))
  {
    double gx = atof(strtok(keyWord+7," "));
    double gy = atof(strtok(nullptr," "));
    double gz = Dim::dimension == 3 ? atof(strtok(nullptr," ")) : 0.0;
    IFEM::cout <<"\nGravitation vector: "<< gx <<" "<< gy;
    if (Dim::dimension == 3) IFEM::cout <<" "<< gz;
    IFEM::cout << std::endl;
    this->getIntegrand()->setGravity(gx,gy,gz);
  }

  else if (!strncasecmp(keyWord,"CONSTANT_PRESSURE",17))
    nConstPress  = atoi(keyWord+17);
  else if (!strncasecmp(keyWord,"LINEAR_PRESSURE",15))
    nLinearPress = atoi(keyWord+15);


  // The remaining keywords are retained for backward compatibility with the
  // prototype version. They enable direct specification of properties onto
  // the topological entities (blocks and faces) of the model.

  else if (!strncasecmp(keyWord,"PRESSURE",8))
  {
    Property press;
    press.pcode = Property::NEUMANN;
    press.ldim = Dim::dimension - 1;

    int npres = atoi(keyWord+8);
    IFEM::cout <<"\nNumber of pressures: "<< npres << std::endl;
    for (int i = 0; i < npres && (cline = utl::readLine(is)); i++)
    {
      press.pindx = 1+i;
      press.patch = atoi(strtok(cline," "));

      int pid = this->getLocalPatchIndex(press.patch);
      if (pid < 0) return false;
      if (pid < 1) continue;

      press.lindx = atoi(strtok(nullptr," "));
      if (press.lindx < 1 || press.lindx > 2*Dim::dimension)
      {
        std::cerr <<" *** SIMElasticity::parse: Invalid face index "
                  << (int)press.lindx << std::endl;
        return false;
      }

      if (Dim::mySol && Dim::mySol->getStressSol())
      {
        IFEM::cout <<"\tTraction on P"<< press.patch
                   << (Dim::dimension==3?" F":" E")
                   << (int)press.lindx << std::endl;
        Dim::myTracs[1+i] = new TractionField(*Dim::mySol->getStressSol());
      }
      else
      {
        int pdir = atoi(strtok(nullptr," "));
        double p = atof(strtok(nullptr," "));
        IFEM::cout <<"\tPressure on P"<< press.patch
                   << (Dim::dimension==3?" F":" E")
                   << (int)press.lindx <<" direction "<< pdir <<": ";
        if ((cline = strtok(nullptr," ")))
        {
          const RealFunc* pf = utl::parseRealFunc(cline,p);
          Dim::myTracs[1+i] = new PressureField(pf,pdir);
        }
        else
        {
          IFEM::cout << p;
          Dim::myTracs[1+i] = new PressureField(p,pdir);
        }
        IFEM::cout << std::endl;
      }

      press.patch = pid;
      Dim::myProps.push_back(press);
    }
  }

  else if (!strncasecmp(keyWord,"MATERIAL",8))
  {
    nmat = atoi(keyWord+8);
    IFEM::cout <<"\nNumber of materials: "<< nmat << std::endl;
    ElasticBase* elInt = this->getIntegrand();
    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      IFEM::cout <<"\tMaterial data: ";
      bool planeStrain = Dim::dimension == 2 ? Elastic::planeStrain : true;
      mVec.push_back(elInt->parseMatProp(cline,planeStrain));

      while ((cline = strtok(nullptr," ")))
        if (!strncasecmp(cline,"ALL",3))
          IFEM::cout <<" (for all patches)"<< std::endl;
        else
        {
          int patch = atoi(cline);
          int pid = this->getLocalPatchIndex(patch);
          if (pid < 0) return false;
          if (pid < 1) continue;

          IFEM::cout <<" (for P"<< patch <<")"<< std::endl;
          Dim::myProps.push_back(Property(Property::MATERIAL,
                                          mVec.size()-1,pid,3));
        }
    }
  }

  else if (!strncasecmp(keyWord,"LOCAL_SYSTEM",12))
  {
    size_t i = 12;
    while (i < strlen(keyWord) && isspace(keyWord[i])) i++;
    static_cast<Elasticity*>(this->getIntegrand())->parseLocalSystem(keyWord+i);
  }

  else if (!strncasecmp(keyWord,"ANASOL",6))
    return this->parseAnaSol(keyWord,is);

  else
    return this->Dim::parse(keyWord,is);

  int npres = nConstPress + nLinearPress;
  if (npres > 0)
  {
    IFEM::cout <<"\nNumber of pressures: "<< npres << std::endl;
    for (int i = 0; i < npres && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      int pdir = atoi(strtok(nullptr," "));
      double p = atof(strtok(nullptr," "));
      IFEM::cout <<"\tPressure code "<< code <<" direction "<< pdir
                 <<": "<< p << std::endl;

      this->setPropertyType(code,Property::NEUMANN);

      if (nLinearPress)
      {
        RealFunc* pfl = new ConstTimeFunc(new LinearFunc(p));
        Dim::myTracs[code] = new PressureField(pfl,pdir);
      }
      else
        Dim::myTracs[code] = new PressureField(p,pdir);
    }
  }

  return true;
}


template<class Dim>
bool SIMElasticity<Dim>::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),"postprocessing"))
  {
    const tinyxml2::XMLElement* child = elem->FirstChildElement();
    for (; child && !plotRgd; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"plot_rigid"))
        plotRgd = true;
      else if (!strcasecmp(child->Value(),"strain"))
        Elasticity::wantStrain = true;
  }

  if (strcasecmp(elem->Value(),myContext.c_str()))
    return this->Dim::parse(elem);

  bool result = true;
  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())

    if (!strcasecmp(child->Value(),"isotropic") ||
        !strcasecmp(child->Value(),"texturematerial"))
    {
      if (!strcasecmp(child->Value(),"texturematerial"))
      {
        Real2DMat domain;
        for (ASMbase* pch : Dim::myModel)
          if (pch->getParameterDomain(domain))
            for (size_t d = 0; d < domain.size(); d++)
              if (domain[d].front() != 0.0 || domain[d].back() != 1.0)
              {
                std::cerr <<" *** Texture material requires unit parametric"
                          <<" domain, "<< char('u'+d) <<"0 = "
                          << domain[d].front() <<", "<< char('u'+d) <<"1 = "
                          << domain[d].back() << std::endl;
                return false;
              }
      }
      IFEM::cout <<"  Parsing <"<< child->Value() <<">"<< std::endl;
      int code = this->parseMaterialSet(child,mVec.size());
      IFEM::cout <<"\tMaterial code "<< code <<":";
      utl::getAttribute(child,"planeStrain",Elastic::planeStrain);
      bool planeStrain = Dim::dimension == 2 ? Elastic::planeStrain : true;
      mVec.push_back(this->getIntegrand()->parseMatProp(child,planeStrain));
    }

    else if (!strcasecmp(child->Value(),"bodyforce"))
    {
      IFEM::cout <<"  Parsing <"<< child->Value() <<">"<< std::endl;
      std::string set, type;
      utl::getAttribute(child,"set",set);
      int code = this->getUniquePropertyCode(set,Dim::dimension==3?123:12);
      if (code == 0) utl::getAttribute(child,"code",code);
      if (child->FirstChild() && code > 0)
      {
        utl::getAttribute(child,"type",type,true);
        IFEM::cout <<"\tBodyforce code "<< code;
        if (!type.empty()) IFEM::cout <<" ("<< type <<")";
        VecFunc* f = utl::parseVecFunc(child->FirstChild()->Value(),type);
        if (f) this->setVecProperty(code,Property::BODYLOAD,f);
        IFEM::cout << std::endl;
      }
    }

    else if (!strcasecmp(child->Value(),"boundaryforce"))
    {
      IFEM::cout <<"  Parsing <"<< child->Value() <<">"<< std::endl;
      std::string set;
      int code = 0;
      if (utl::getAttribute(child,"set",set))
        code = this->getUniquePropertyCode(set);
      else if (!utl::getAttribute(child,"code",code) || code == 0)
        continue;

      IFEM::cout <<"\tBoundary force ";
      if (!set.empty()) IFEM::cout <<"\""<< set <<"\" ";
      IFEM::cout <<"code "<< code << std::endl;
      this->setPropertyType(code,Property::OTHER);

      bCode[code] = Vec3();
    }

    else if (!strcasecmp(child->Value(),"rigid"))
      result &= this->parseRigid(child,this);

    else if (!strcasecmp(child->Value(),"anasol"))
      result &= this->parseAnaSol(child);

    else if (!strcasecmp(child->Value(),"dualfield"))
    {
      FunctionBase* exf = this->parseDualTag(child,2);
      if (exf)
        static_cast<Elasticity*>(this->getIntegrand())->addExtrFunction(exf);
    }
    else if (!this->getIntegrand()->parse(child))
      result &= this->Dim::parse(child);

  return result;
}


/*!
  This method is reimplemented to handle dirichlet conditions on the explicit
  master nodes of rigid couplings which not are regular nodes in a patch.
  These nodes may also have rotational degrees of freedom.
*/

template<class Dim>
bool SIMElasticity<Dim>::addConstraint (int patch, int lndx, int ldim,
                                        int dirs, int code, int& ngnod,
                                        char basis)
{
  if (patch == 0 && ldim == 0)
  {
    typename Dim::IdxVec3* masterPt = this->getDiscretePoint(lndx);
    if (masterPt) {
      for (ASMbase* pch : Dim::myModel)
      {
        // Check if this patch has master points that should be constrained.
        // Must increase the number of field variables temporarily to account
        // for the rotational degrees of freedom
        unsigned char oldnf = pch->getNoFields();
        pch->setNoFields(pch->getNoSpaceDim()*(pch->getNoSpaceDim()+1)/2);
        bool found = pch->constrainXnode(masterPt->first,dirs,code);
        pch->setNoFields(oldnf);
        if (found) return true;
      }

      std::cerr <<" *** SIMElasticity::addConstraint: Master point "
                << masterPt->first <<" ("<< masterPt->second
                <<") is not present in any patch.\n";
      return false;
    }
  }

  return this->Dim::addConstraint(patch,lndx,ldim,dirs,code,ngnod,basis);
}


template<class Dim>
bool SIMElasticity<Dim>::initMaterial (size_t propInd)
{
  Elasticity* elp = dynamic_cast<Elasticity*>(Dim::myProblem);
  if (!elp) return false;

  if (propInd >= mVec.size()) propInd = mVec.size()-1;

  elp->setMaterial(mVec[propInd]);
  return true;
}


template<class Dim>
bool SIMElasticity<Dim>::initBodyLoad (size_t patchInd)
{
  Elasticity* elp = dynamic_cast<Elasticity*>(Dim::myProblem);
  if (!elp) return false;

  elp->setBodyForce(this->getVecFunc(patchInd,Property::BODYLOAD));
  return true;
}


template<class Dim>
bool SIMElasticity<Dim>::initNeumann (size_t propInd)
{
  Elasticity* elp = dynamic_cast<Elasticity*>(Dim::myProblem);
  if (!elp) return false;

  typename Dim::VecFuncMap::const_iterator vit = Dim::myVectors.find(propInd);
  typename Dim::TracFuncMap::const_iterator tit = Dim::myTracs.find(propInd);

  if (vit != Dim::myVectors.end())
    elp->setTraction(vit->second);
  else if (tit != Dim::myTracs.end())
    elp->setTraction(tit->second);
  else
    return false;

  return true;
}


template<class Dim>
size_t SIMElasticity<Dim>::getVolumeIndex () const
{
  Elasticity* elp = dynamic_cast<Elasticity*>(Dim::myProblem);
  return elp ? elp->getVolumeIndex(this->haveAnaSol()) : 0;
}


template<class Dim>
bool SIMElasticity<Dim>::postProcessNorms (Vectors& gNorm, Matrix* eNorm)
{
  return this->revertSqrt(gNorm,eNorm);
}


template<class Dim>
void SIMElasticity<Dim>::printNormGroup (const Vector& gNorm,
                                         const Vector& rNorm,
                                         const std::string& prjName) const
{
  Elastic::printNorms(gNorm,rNorm,prjName,this);
}


/*!
  This method is reimplemented to account for potential rigid couplings.
*/

template<class Dim>
bool SIMElasticity<Dim>::writeGlvG (int& nBlock,
                                    const char* inpFile, bool doClear)
{
  if (!this->Dim::writeGlvG(nBlock,inpFile,doClear))
    return false;
  else if (!plotRgd)
    return true;

  ElementBlock* rgd = this->rigidGeometry(this);
  return rgd ? this->getVTF()->writeGrid(rgd,"Rigid couplings",++nBlock) : true;
}


/*!
  The boundaries for which the interface forces are extracted and printed
  are identified by the property set codes in \ref bCode, which are assigned
  assigned values by parsing the `<boundaryforce>` tags in the input file.
*/

template<class Dim>
void SIMElasticity<Dim>::printIFforces (const Vector& sf, RealArray& weights)
{
  Elastic::printBoundaryForces(sf,weights,bCode,this);
}


template class SIMElasticity<SIM2D>;
template class SIMElasticity<SIM3D>;
