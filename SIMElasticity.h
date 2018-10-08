// $Id$
//==============================================================================
//!
//! \file SIMElasticity.h
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#ifndef _SIM_ELASTICITY_H
#define _SIM_ELASTICITY_H

#include "IFEM.h"
#include "ASMs2D.h"
#include <GoTools/geometry/SplineSurface.h>
#include "Elasticity.h"
#include "ElasticityUtils.h"
#include "MaterialBase.h"
#include "IsotropicTextureMat.h"
#include "ForceIntegrator.h"
#include "Property.h"
#include "TimeStep.h"
#include "AnaSol.h"
#include "Functions.h"
#include "Utilities.h"
#include "tinyxml.h"

typedef std::vector<Material*> MaterialVec; //!< Convenience declaration


/*!
  \brief Driver class for isogeometric FEM analysis of elasticity problems.
  \details The class incapsulates data and methods for solving elasticity
  problems using NURBS-based finite elements. It reimplements the parse methods
  and some property initialization methods of the parent class.
*/

template<class Dim> class SIMElasticity : public Dim
{
public:
  //! \brief Default constructor.
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  explicit SIMElasticity(bool checkRHS = false) : Dim(Dim::dimension,checkRHS)
  {
    myContext = "elasticity";
    aCode = bCode = 0;
  }

  //! \brief The destructor frees the dynamically allocated material properties.
  virtual ~SIMElasticity()
  {
    // To prevent the SIMbase destructor try to delete already deleted functions
    if (aCode > 0)
      Dim::myVectors.erase(aCode);

    for (Material* mat : mVec)
      delete mat;
  }

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const { return "Elasticity"; }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp)
  {
    Elasticity* elp = dynamic_cast<Elasticity*>(Dim::myProblem);
    if (elp)
      elp->advanceStep(tp.time.dt,tp.time.dtn);

    return true;
  }

  //! \brief Initializes the property containers of the model.
  virtual void clearProperties()
  {
    // To prevent SIMbase::clearProperties deleting the analytical solution
    if (aCode > 0)
      Dim::myVectors.erase(aCode);
    aCode = 0;

    Elasticity* elp = dynamic_cast<Elasticity*>(Dim::myProblem);
    if (elp)
    {
      elp->setMaterial(nullptr);
      elp->setBodyForce(nullptr);
      elp->setTraction((VecFunc*)nullptr);
      elp->setTraction((TractionFunc*)nullptr);
    }

    for (Material* mat : mVec)
      delete mat;
    mVec.clear();

    this->Dim::clearProperties();
  }

  //! \brief Calculates the traction resultant associated with a given boundary.
  //! \param[out] f Calculated traction resultant
  //! \param[in] sol Primary solution vectors
  //! \param[in] tp Time stepping parameters
  //!
  //! \details The boundary for which the traction is calculated is identified
  //! by the property set code \a bCode which is assigned value by parsing
  //! the `<boundaryforce>` tag in the input file.
  bool getBoundaryForce(Vector& f, const Vectors& sol, const TimeStep& tp)
  {
    if (bCode == 0) return false;

    f = SIM::getBoundaryForce(sol,this,bCode,tp.time);
    Dim::adm.allReduceAsSum(f);
    return true;
  }

  //! \brief Extracts the reaction forces associated with a given boundary.
  //! \param[out] rf Reaction force resultant for specified boundary
  //!
  //! \details The boundary for which the reaction force is returned
  //! is identified by the property set code \a bCode which is assigned value
  //! by parsing the `<boundaryforce>` tag in the input file.
  bool getBoundaryReactions(Vector& rf)
  {
    if (bCode == 0) return false;

    bool ok = this->getCurrentReactions(rf,bCode);
    Dim::adm.allReduceAsSum(rf);
    return ok;
  }

  //! \brief Returns whether an analytical solution is available or not.
  virtual bool haveAnaSol() const
  {
    return (Dim::mySol && Dim::mySol->getStressSol());
  }

protected:
  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented inserting a call to \a getIntegrand.
  //! This makes sure the integrand has been allocated in case of minimum input.
  //! It also resolves inhomogeneous boundary condition fields in case they are
  //! derived from the analytical solution.
  virtual void preprocessA()
  {
    this->getIntegrand();
    this->printProblem();

    // Deactivate principal stress ouput for Lagrange/Spectral interpolations
    if (Dim::opt.discretization < ASM::Spline)
      Elasticity::wantPrincipalStress = false;

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
          Dim::myTracs[p.pindx] = new TractionField(*stressField);
        }
        else
          p.pcode = Property::UNDEFINED;
      }
  }

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to ensure that threading groups are
  //! established for the patch faces subjected to boundary force integration.
  virtual bool preprocessB()
  {
    if (bCode == 0) return true;

    for (const Property& p : Dim::myProps)
      if (bCode == p.pindx)
        this->generateThreadGroups(p,Dim::msgLevel < 2);

    return true;
  }

  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand() = 0;
  //! \brief Parses a dimension-specific data section from an input file.
  virtual bool parseDimSpecific(char*, std::istream&) { return false; }
  //! \brief Parses a dimension-specific data section from an XML element.
  virtual bool parseDimSpecific(const TiXmlElement*) { return false; }

  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is)
  {
    char* cline = nullptr;
    int nmat = 0;
    int nConstPress = 0;
    int nLinearPress = 0;

    if (this->parseDimSpecific(keyWord,is))
      return true;

    else if (!strncasecmp(keyWord,"ISOTROPIC",9))
    {
      nmat = atoi(keyWord+10);
      IFEM::cout <<"\nNumber of isotropic materials: "<< nmat << std::endl;
      Elasticity* elInt = this->getIntegrand();
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
      IFEM::cout <<"\nGravitation vector: " << gx <<" "<< gy;
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
          std::cerr <<" *** SIMElasticity3D::parse: Invalid face index "
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
      Elasticity* elInt = this->getIntegrand();
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
      this->getIntegrand()->parseLocalSystem(keyWord+i);
    }

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

  //! \brief Parses a data section from an XML element
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),myContext.c_str()))
      return this->Dim::parse(elem);

    bool result = true;
    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (this->parseDimSpecific(child))
        continue;
      else if (!strcasecmp(child->Value(),"texturematerial"))
      {
        if (Dim::dimension != 2)
        {
          std::cerr << "Texture material not supported for trivariate models" << std::endl;
          return false;
        }
        for (size_t i=1; i<=this->getNoPatches(); ++i)
        {
          ASMs2D *patch = dynamic_cast<ASMs2D*>(this->getPatch(i));
          if (!patch)
          {
            std::cerr << "Only works for ASMs2D..." << std::endl;
            return false;
          }
          Go::SplineSurface *surf = patch->getSurface();
          if (surf->startparam_u() != 0 || surf->endparam_u() != 1 ||
              surf->startparam_v() != 0 || surf->endparam_v() != 1)
          {
            std::cerr << "Texture material requires unit parametric domain" << std::endl;
            return false;
          }
        }
        IFEM::cout <<"  Parsing <"<< child->Value() <<">"<< std::endl;
        int code = this->parseMaterialSet(child,mVec.size());
        IFEM::cout <<"\tMaterial code "<< code <<":";
        bool planeStrain = Dim::dimension == 2 ? Elastic::planeStrain : true;
        IsotropicTextureMat* mat = new IsotropicTextureMat(planeStrain,
                                                           Elastic::axiSymmetry);
        mat->parse(child);
        mVec.push_back(mat);
      }
      else if (!strcasecmp(child->Value(),"isotropic"))
      {
        IFEM::cout <<"  Parsing <"<< child->Value() <<">"<< std::endl;
        int code = this->parseMaterialSet(child,mVec.size());
        IFEM::cout <<"\tMaterial code "<< code <<":";
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
        if (utl::getAttribute(child,"set",set))
          bCode = this->getUniquePropertyCode(set);
        else if (!utl::getAttribute(child,"code",bCode) || bCode == 0)
          continue;

        IFEM::cout <<"\tBoundary force ";
        if (!set.empty()) IFEM::cout <<"\""<< set <<"\" ";
        IFEM::cout <<"code "<< bCode << std::endl;
        this->setPropertyType(bCode,Property::OTHER);
      }

      else if (!this->getIntegrand()->parse(child))
        result &= this->Dim::parse(child);

    return result;
  }

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd)
  {
    Elasticity* elp = dynamic_cast<Elasticity*>(Dim::myProblem);
    if (!elp) return false;

    if (propInd >= mVec.size()) propInd = mVec.size()-1;

    elp->setMaterial(mVec[propInd]);
    return true;
  }

  //! \brief Initializes the body load properties for current patch.
  //! \param[in] patchInd 1-based patch index
  virtual bool initBodyLoad(size_t patchInd)
  {
    Elasticity* elp = dynamic_cast<Elasticity*>(Dim::myProblem);
    if (!elp) return false;

    elp->setBodyForce(this->getVecFunc(patchInd,Property::BODYLOAD));
    return true;
  }

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd)
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

public:
  //! \brief Prints a norm group to the log stream.
  //! \param[in] gNorm The norm values to print
  //! \param[in] rNorm Reference norms for the first norm group
  //! \param[in] prjName Projection name associated with this norm group
  virtual void printNormGroup (const Vector& gNorm, const Vector& rNorm,
                               const std::string& prjName) const
  {
    Elastic::printNorms(gNorm,rNorm,prjName,this);
  }

protected:
  MaterialVec mVec;      //!< Material data
  std::string myContext; //!< XML-tag to search for problem inputs within

private:
  int aCode; //!< Analytical BC code (used by destructor)
  int bCode; //!< Property set code for boundary traction resultant calculation
};

#endif
