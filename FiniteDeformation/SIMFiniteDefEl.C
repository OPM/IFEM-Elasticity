// $Id$
//==============================================================================
//!
//! \file SIMFiniteDefEl.C
//!
//! \date Aug 22 2025
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution drivers for NURBS-based finite deformation analysis.
//!
//==============================================================================

#include "SIMFiniteDefEl.h"
#include "LinIsotropic.h"
#include "LinearMaterial.h"
#include "NeoHookeMaterial.h"
#include "PlasticMaterial.h"
#include "DruckerPrager.h"
#include "ElasticBase.h"
#include "ElasticityUtils.h"

#include "IFEM.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "AlgEqSystem.h"
#include "ASMmxBase.h"
#include "Functions.h"
#include "Utilities.h"
#include "Property.h"
#include "tinyxml2.h"


template<class Dim>
SIMFiniteDefEl<Dim>::SIMFiniteDefEl (bool rhs, const std::vector<int>& options)
  : SIMElasticity<Dim>(rhs), mDat(SIMElasticity<Dim>::mVec), nlo(Dim::dimension)
{
  nlo.form = options.size() > 0 ? options[0] : SIM::TOTAL_LAGRANGE;
  nlo.pOrd = options.size() > 1 ? options[1] : 0;
  if (nlo.form == SIM::MIXED_QnQn1)
    Dim::nf.push_back(2);
}


template<class Dim>
bool SIMFiniteDefEl<Dim>::parse (char* keyWord, std::istream& is)
{
  ElasticBase* elInt = this->getIntegrand();
  Material* defaultMat = nullptr;

  if (!strncasecmp(keyWord,"ISOTROPIC",9))
  {
    int nmat = atoi(keyWord+9);
    IFEM::cout <<"\nNumber of isotropic materials: "<< nmat << std::endl;

    char* cline = nullptr;
    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      if (code > 0)
        this->setPropertyType(code,Property::MATERIAL,mDat.size());

      double E   = atof(strtok(nullptr," "));
      double nu  = atof(strtok(nullptr," "));
      double rho = atof(strtok(nullptr," "));
      int matVer = (cline = strtok(nullptr," ")) ? atoi(cline) : -1;
      if (matVer >= 0 && nlo.form >= SIM::UPDATED_LAGRANGE)
        mDat.push_back(new NeoHookeMaterial(E,nu,rho,matVer));
      else if (Dim::dimension == 2)
        mDat.push_back(new LinIsotropic(E,nu,rho,
                                        !Elastic::planeStrain,
                                        Elastic::axiSymmetry));
      else
        mDat.push_back(new LinIsotropic(E,nu,rho));

      if (matVer < 0 && nlo.form >= SIM::UPDATED_LAGRANGE)
        mDat.back() = new LinearMaterial(mDat.back());

      IFEM::cout <<"\tMaterial code "<< code <<": "<< E <<" "<< nu
                 <<" "<< rho <<" ("<< matVer <<")"<< std::endl;
      if (!defaultMat || code == 0) defaultMat = mDat.back();
    }
  }

  else if (!strncasecmp(keyWord,"PLASTIC",7))
  {
    int nmat = atoi(keyWord+7);
    IFEM::cout <<"\nNumber of plastic materials: "<< nmat << std::endl;

    char* cline = nullptr;
    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      if (code > 0)
        this->setPropertyType(code,Property::MATERIAL,mDat.size());

      RealArray pMAT;
      while ((cline = strtok(nullptr," ")))
        pMAT.push_back(atof(cline));
      mDat.push_back(new PlasticMaterial(pMAT));

      IFEM::cout <<"\tMaterial code "<< code <<":";
      for (double v : pMAT) IFEM::cout <<" "<< v;
      IFEM::cout << std::endl;
      if (!defaultMat || code == 0) defaultMat = mDat.back();
    }
  }

  else if (!strncasecmp(keyWord,"MATERIAL",8))
  {
    std::cerr <<" *** SIMFiniteDefEl::parse: The keyword MATERIAL"
              <<" is not supported\n     for finite deformation analysis."
              <<" You must use ISOTROPIC instead."<< std::endl;
    return false;
  }

  else
    return this->SIMElasticity<Dim>::parse(keyWord,is);

  if (elInt && defaultMat)
    elInt->setMaterial(defaultMat);

  return true;
}


template<class Dim>
bool SIMFiniteDefEl<Dim>::parse (const tinyxml2::XMLElement* elem)
{
  if (strcasecmp(elem->Value(),"finitedeformation"))
    return this->SIMElasticity<Dim>::parse(elem);

  ElasticBase* elInt = this->getIntegrand();
  Material* defaultMat = nullptr;

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())

    if (!strcasecmp(child->Value(),"isotropic"))
    {
      int code = this->parseMaterialSet(child,mDat.size());
      IFEM::cout <<"\tMaterial code "<< code;

      int mVer = -1;
      if (utl::getAttribute(child,"version",mVer))
        IFEM::cout <<" ("<< mVer <<"):";
      else
        IFEM::cout <<":";

      if (mVer >= 0 && nlo.form >= SIM::UPDATED_LAGRANGE)
        mDat.push_back(new NeoHookeMaterial(mVer));
      else if (Dim::dimension == 2)
        mDat.push_back(new LinIsotropic(!Elastic::planeStrain,
                                        Elastic::axiSymmetry));
      else
        mDat.push_back(new LinIsotropic());

      mDat.back()->parse(child);
      IFEM::cout << std::endl;

      if (mVer < 0 && nlo.form >= SIM::UPDATED_LAGRANGE)
        mDat.back() = new LinearMaterial(mDat.back());

      if (!defaultMat || code == 0) defaultMat = mDat.back();
    }

    else if (!strcasecmp(child->Value(),"druckerprager"))
    {
      int code = this->parseMaterialSet(child,mDat.size());
      IFEM::cout <<"\tMaterial code "<< code;

      mDat.push_back(new DruckerPrager(Dim::dimension,
                                       !Elastic::planeStrain,
                                       Elastic::axiSymmetry));
      mDat.back()->parse(child);
      IFEM::cout << std::endl;

      if (nlo.form >= SIM::UPDATED_LAGRANGE)
        mDat.back() = new LinearMaterial(mDat.back());

      if (!defaultMat || code == 0) defaultMat = mDat.back();
    }
    else if (!strcasecmp(child->Value(),"plastic"))
    {
      int code = this->parseMaterialSet(child,mDat.size());

      RealArray pMAT;
      ScalarFunc* hfunc = nullptr;
      const tinyxml2::XMLElement* gchild = child->FirstChildElement();
      for (; gchild; gchild = gchild->NextSiblingElement())
        if (strcasecmp(gchild->Value(),"hardeningcurve"))
        {
          std::string value(gchild->Value());
          char* cval = strtok(const_cast<char*>(value.c_str())," ");
          for (; cval; cval = strtok(nullptr," "))
            pMAT.push_back(atof(cval));
        }
        else if (gchild->FirstChild())
        {
          // An isotropic hardening function is specified
          IFEM::cout <<" Hardening: ";
          std::string type("expression");
          utl::getAttribute(gchild,"type",type);
          hfunc = utl::parseTimeFunc(gchild->FirstChild()->Value(),type);
        }

      if (pMAT.size() < 11) pMAT.resize(11,0.0);
      utl::getAttribute(child,"Bmod" ,pMAT[0]);
      utl::getAttribute(child,"Emod" ,pMAT[0]);
      utl::getAttribute(child,"E"    ,pMAT[0]);
      utl::getAttribute(child,"Smod" ,pMAT[1]);
      utl::getAttribute(child,"nu"   ,pMAT[1]);
      utl::getAttribute(child,"rho"  ,pMAT[3]);
      utl::getAttribute(child,"Hiso" ,pMAT[4]);
      utl::getAttribute(child,"Hkin" ,pMAT[5]);
      utl::getAttribute(child,"yield",pMAT[6]);
      utl::getAttribute(child,"Y0"   ,pMAT[7]);
      utl::getAttribute(child,"Yinf" ,pMAT[8]);
      utl::getAttribute(child,"beta" ,pMAT[9]);
      utl::getAttribute(child,"istrt",pMAT[10]);
      int iYield = static_cast<int>(pMAT[6]);
      if (iYield == 4)
      {
        utl::getAttribute(child,"A"  ,pMAT[7]);
        utl::getAttribute(child,"B"  ,pMAT[8]);
        utl::getAttribute(child,"n"  ,pMAT[9]);
      }
      else if (iYield == 6)
      {
        if (pMAT.size() < 13) pMAT.resize(13,0.0);
        utl::getAttribute(child,"Q1",pMAT[8]);
        utl::getAttribute(child,"C1",pMAT[9]);
        utl::getAttribute(child,"Q2",pMAT[11]);
        utl::getAttribute(child,"C2",pMAT[12]);
      }
      mDat.push_back(new PlasticMaterial(pMAT,hfunc));

      IFEM::cout <<"\tMaterial code "<< code <<":";
      for (double v : pMAT) IFEM::cout <<" "<< v;
      IFEM::cout << std::endl;

      if (!defaultMat || code == 0) defaultMat = mDat.back();
    }

    else if (!strcasecmp(child->Value(),"contact"))
    {
      if (!this->createFEMmodel(false))
        return false;

      if (!this->parseContactTag(child,Dim::myModel,Dim::myEntitys))
        return false;

      Dim::preserveNOrder = true; // because extra nodes have been added
    }

    else if (elInt)
      elInt->parse(child);

  if (elInt && defaultMat)
    elInt->setMaterial(defaultMat);

  return true;
}


template<class Dim>
ElasticBase* SIMFiniteDefEl<Dim>::getIntegrand ()
{
  if (!Dim::myProblem) Dim::myProblem = nlo.getIntegrand();
  return dynamic_cast<ElasticBase*>(Dim::myProblem);
}


template<class Dim>
bool SIMFiniteDefEl<Dim>::preprocessB ()
{
  Dim::myInts.insert({ 0, Dim::myProblem });
  if (this->withContact())
  {
    Dim::opt.num_threads_SLU *= -1; // do not lock the sparsity pattern
    return this->preprocessContact(Dim::myInts,*this->getSAM(),
                                   this->getNoSpaceDim());
  }

  return true;
}


template<class Dim>
bool SIMFiniteDefEl<Dim>::createContactSet (const std::string& slaveSet,
                                            int& code)
{
  if (!(code = this->getUniquePropertyCode(slaveSet)))
    return false;

  this->setPropertyType(code,Property::NEUMANN_GENERIC);
  return true;
}


template<class Dim>
bool SIMFiniteDefEl<Dim>::preprocessBeforeAsmInit (int& ngnod)
{
  this->renumberContactBodies(*Dim::g2l);
  for (ASMbase* pch : Dim::myModel)
    this->addLagrangeMultipliers(pch,ngnod);

  return this->SIMElasticity<Dim>::preprocessBeforeAsmInit(ngnod);
}


template<class Dim>
bool SIMFiniteDefEl<Dim>::assembleDiscreteTerms (const IntegrandBase* problem,
                                                 const TimeDomain&)
{
  if (!Dim::myEqSys) return true;

  return this->assembleMortarTangent(problem,
                                     Dim::myEqSys->getMatrix(),
                                     Dim::myEqSys->getVector());
}


template<class Dim>
bool SIMFiniteDefEl<Dim>::printProblem () const
{
  if (!this->SIMElasticity<Dim>::printProblem())
    return false;
  else if (Dim::opt.discretization < ASM::Spline)
    return true;

  if (ASMmxBase::Type == ASMmxBase::FULL_CONT_RAISE_BASIS1 ||
      ASMmxBase::Type == ASMmxBase::FULL_CONT_RAISE_BASIS2)
    IFEM::cout <<"Using C^(p-1) continuous displacement basis\n";
  else if (nlo.form == SIM::MIXED_QnQn1)
    IFEM::cout <<"Using C^(p-2) continuous displacement basis\n";

  return true;
}


template<class Dim>
bool SIMFiniteDefEl<Dim>::updateDirichlet (double time, const Vector* prevSol)
{
  if (prevSol)
    this->SIMContact::updateDirichlet(time);

  return this->SIMElasticity<Dim>::updateDirichlet(time,prevSol);
}


template<class Dim>
bool SIMFiniteDefEl<Dim>::updateConfiguration (const Vector& solution)
{
  return this->updateContactBodies(solution);
}


template<class Dim>
void SIMFiniteDefEl<Dim>::printIFforces (const Vector& sf, RealArray& weights)
{
  Elastic::printBoundaryForces(sf,weights,this->bCode,this,true);
}


template<class Dim>
bool SIMFiniteDefEl<Dim>::writeGlvG (int& nBlock, const char* inpFile, bool)
{
  if (!this->SIMElasticity<Dim>::writeGlvG(nBlock,inpFile))
    return false;

  return this->writeGlvBodies(this->getVTF(),nBlock);
}


template<class Dim>
bool SIMFiniteDefEl<Dim>::writeGlvA (int& nBlock, int iStep, double, int) const
{
  return this->writeGlvBodyMovements(this->getVTF(),iStep,nBlock);
}


template<class Dim>
void SIMFiniteDefEl<Dim>::dumpMoreResults (double, utl::LogStream& os,
                                           std::streamsize prec) const
{
  const RealArray* rf = Dim::myEqSys ? Dim::myEqSys->getReactions() : nullptr;
  if (rf) this->printBodyReactions(*this->getSAM(),*rf,os,prec);
}


template class SIMFiniteDefEl<SIM2D>;
template class SIMFiniteDefEl<SIM3D>;
