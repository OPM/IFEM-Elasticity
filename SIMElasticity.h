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

#include "SIM2D.h"
#include "SIM3D.h"
#include "Elasticity.h"
#include "LinIsotropic.h"
#include "AnaSol.h"
#include "Functions.h"
#include "Utilities.h"
#include "tinyxml.h"
#include "TimeStep.h"

/*!
  \brief Driver class for isogeometric FEM analysis of elasticity problems.
  \details The class incapsulates data and methods for solving linear elasticity
  problems using NURBS-based finite elements. It reimplements the parse methods
  and some property initialization methods of the parent class.
*/

template<class Dim> class SIMElasticity : public Dim
{
public:
  //! \brief Default constructor.
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  SIMElasticity(bool checkRHS = false) : Dim(Dim::dimension,0,checkRHS)
  {
    aCode = 0;
  }

  //! \brief The destructor frees the dynamically allocated material properties.
  virtual ~SIMElasticity()
  {
    // To prevent the SIMbase destructor try to delete already deleted functions
    if (aCode > 0)
      Dim::myVectors.erase(aCode);

    for (size_t i = 0; i < mVec.size(); i++)
      delete mVec[i];
  }

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
      elp->setMaterial(NULL);
      elp->setBodyForce(NULL);
      elp->setTraction((VecFunc*)NULL);
      elp->setTraction((TractionFunc*)NULL);
    }

    for (size_t i = 0; i < mVec.size(); i++)
      delete mVec[i];
    mVec.clear();

    this->Dim::clearProperties();
    this->Dim::setOwnProblem(true); // we always own our integrand
  }

protected:
  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented inserting a call to \a getIntegrand.
  //! This makes sure the integrand has been allocated in case of minimum input.
  //! It also resolves inhomogeneous boundary condition fields in case they are
  //! derived from the analytical solution.
  virtual void preprocessA()
  {
    if (!Dim::myProblem)
    {
      this->getIntegrand();
      if (Dim::myPid == 0)
        this->printProblem(std::cout);
    }

    if (!Dim::mySol) return;

    // Define analytical boundary condition fields
    PropertyVec::iterator p;
    for (p = Dim::myProps.begin(); p != Dim::myProps.end(); ++p)
      if (p->pcode == Property::DIRICHLET_ANASOL)
      {
        VecFunc* vecField = Dim::mySol->getVectorSol();
        if (!vecField)
          p->pcode = Property::UNDEFINED;
        else if (aCode == abs(p->pindx))
          p->pcode = Property::DIRICHLET_INHOM;
        else if (aCode == 0)
        {
          aCode = abs(p->pindx);
          Dim::myVectors[aCode] = vecField;
          p->pcode = Property::DIRICHLET_INHOM;
        }
        else
          p->pcode = Property::UNDEFINED;
      }
      else if (p->pcode == Property::NEUMANN_ANASOL)
      {
        STensorFunc* stressField = Dim::mySol->getStressSol();
        if (stressField)
        {
          p->pcode = Property::NEUMANN;
          Dim::myTracs[p->pindx] = new TractionField(*stressField);
        }
        else
          p->pcode = Property::UNDEFINED;
      }
  }

public:
  static bool planeStrain; //!< Plane strain/stress option - 2D only
  static bool axiSymmetry; //!< Axisymmtry option - 2D only

protected:
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand() { return NULL; }

  //! \brief Parses a dimension-specific data section from an input file.
  virtual bool parseDimSpecific(char*, std::istream&) { return false; }
  //! \brief Parses a dimension-specific data section from an XML element.
  virtual bool parseDimSpecific(const TiXmlElement*) { return false; }

  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is)
  {
    char* cline = NULL;
    int nmat = 0;
    int nConstPress = 0;
    int nLinearPress = 0;

    if (this->parseDimSpecific(keyWord,is))
      return true;

    else if (!strncasecmp(keyWord,"ISOTROPIC",9))
    {
      nmat = atoi(keyWord+10);
      std::cout <<"\nNumber of isotropic materials: "<< nmat << std::endl;

      for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
      {
        int code = atoi(strtok(cline," "));
        if (code > 0)
          this->setPropertyType(code,Property::MATERIAL,mVec.size());

        double E   = atof(strtok(NULL," "));
        double nu  = atof(strtok(NULL," "));
        double rho = atof(strtok(NULL," "));
        if (Dim::dimension == 2)
          mVec.push_back(new LinIsotropic(E,nu,rho,!planeStrain,axiSymmetry));
        else
          mVec.push_back(new LinIsotropic(E,nu,rho));
        std::cout <<"\tMaterial code "<< code <<": "
                  << E <<" "<< nu <<" "<< rho << std::endl;
      }
    }

    else if (!strncasecmp(keyWord,"GRAVITY",7))
    {
      double gx = atof(strtok(keyWord+7," "));
      double gy = atof(strtok(NULL," "));
      double gz = Dim::dimension == 3 ? atof(strtok(NULL," ")) : 0.0;
      std::cout <<"\nGravitation vector: " << gx <<" "<< gy;
      if (Dim::dimension == 3) std::cout <<" "<< gz;
      std::cout << std::endl;
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
      std::cout <<"\nNumber of pressures: "<< npres << std::endl;
      for (int i = 0; i < npres && (cline = utl::readLine(is)); i++)
      {
        press.pindx = 1+i;
        press.patch = atoi(strtok(cline," "));

        int pid = this->getLocalPatchIndex(press.patch);
        if (pid < 0) return false;
        if (pid < 1) continue;

        press.lindx = atoi(strtok(NULL," "));
        if (press.lindx < 1 || press.lindx > 2*Dim::dimension)
        {
          std::cerr <<" *** SIMElasticity3D::parse: Invalid face index "
                    << (int)press.lindx << std::endl;
          return false;
        }

        if (Dim::mySol && Dim::mySol->getStressSol())
        {
          std::cout <<"\tTraction on P"<< press.patch
                    << (Dim::dimension==3?" F":" E")
                    << (int)press.lindx << std::endl;
          Dim::myTracs[1+i] = new TractionField(*Dim::mySol->getStressSol());
        }
        else
        {
          int pdir = atoi(strtok(NULL," "));
          double p = atof(strtok(NULL," "));
          std::cout <<"\tPressure on P"<< press.patch
                    << (Dim::dimension==3?" F":" E")
                    << (int)press.lindx <<" direction "<< pdir <<": ";
          if ((cline = strtok(NULL," ")))
          {
            const RealFunc* pf = utl::parseRealFunc(cline,p);
            Dim::myTracs[1+i] = new PressureField(pf,pdir);
          }
          else
          {
            std::cout << p;
            Dim::myTracs[1+i] = new PressureField(p,pdir);
          }
          std::cout << std::endl;
        }

        press.patch = pid;
        Dim::myProps.push_back(press);
      }
    }

    else if (!strncasecmp(keyWord,"MATERIAL",8))
    {
      nmat = atoi(keyWord+8);
      std::cout <<"\nNumber of materials: "<< nmat << std::endl;
      for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
      {
        double E   = atof(strtok(cline," "));
        double nu  = atof(strtok(NULL," "));
        double rho = atof(strtok(NULL," "));
        while ((cline = strtok(NULL," ")))
          if (!strncasecmp(cline,"ALL",3))
          {
            std::cout <<"\tMaterial for all patches: "
                      << E <<" "<< nu <<" "<< rho << std::endl;
            if (Dim::dimension == 2)
              mVec.push_back(new LinIsotropic(E,nu,rho,
                                              !planeStrain,axiSymmetry));
            else
              mVec.push_back(new LinIsotropic(E,nu,rho));
          }
          else
          {
            int patch = atoi(cline);
            int pid = this->getLocalPatchIndex(patch);
            if (pid < 0) return false;
            if (pid < 1) continue;

            std::cout <<"\tMaterial for P"<< patch
                      <<": "<< E <<" "<< nu <<" "<< rho << std::endl;
            Dim::myProps.push_back(Property(Property::MATERIAL,
                                            mVec.size(),pid,3));
            if (Dim::dimension == 2)
              mVec.push_back(new LinIsotropic(E,nu,rho,
                                              !planeStrain,axiSymmetry));
            else
              mVec.push_back(new LinIsotropic(E,nu,rho));
          }
      }
    }

    else
      return this->Dim::parse(keyWord,is);

    int npres = nConstPress + nLinearPress;
    if (npres > 0)
    {
      std::cout <<"\nNumber of pressures: "<< npres << std::endl;
      for (int i = 0; i < npres && (cline = utl::readLine(is)); i++)
      {
        int code = atoi(strtok(cline," "));
        int pdir = atoi(strtok(NULL," "));
        double p = atof(strtok(NULL," "));
        std::cout <<"\tPressure code "<< code <<" direction "<< pdir
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

    if (nmat > 0 && !mVec.empty())
      this->getIntegrand()->setMaterial(mVec.front());

    return true;
  }

  //! \brief Parses a data section from an XML element
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"elasticity"))
      return this->Dim::parse(elem);

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement()) {
      if (this->parseDimSpecific(child))
        continue;

      if (!strcasecmp(child->Value(),"isotropic")) {
        int code = this->parseMaterialSet(child,mVec.size());

        double E = 1000.0, nu = 0.3, rho = 1.0;
        utl::getAttribute(child,"E",E);
        utl::getAttribute(child,"nu",nu);
        utl::getAttribute(child,"rho",rho);

        if (Dim::dimension == 2)
          mVec.push_back(new LinIsotropic(E,nu,rho,!planeStrain,axiSymmetry));
        else
          mVec.push_back(new LinIsotropic(E,nu,rho));

        std::cout <<"\tMaterial code "<< code <<": "
                  << E <<" "<< nu <<" "<< rho << std::endl;
      }

      else if (!strcasecmp(child->Value(),"gravity")) {
        double gx = 0.0, gy = 0.0, gz = 0.0;
        utl::getAttribute(child,"x",gx);
        utl::getAttribute(child,"y",gy);
        if (Dim::dimension == 3)
          utl::getAttribute(child,"z",gz);
        this->getIntegrand()->setGravity(gx,gy,gz);

        std::cout <<"\tGravitation vector: "<< gx <<" "<< gy;
        if (Dim::dimension == 3) std::cout <<" "<< gz;
        std::cout << std::endl;
      }

      else if (!strcasecmp(child->Value(),"bodyforce")) {
        std::string set, type;
        utl::getAttribute(child,"set",set);
        int code = this->getUniquePropertyCode(set,Dim::dimension==3?123:12);
        if (code == 0) utl::getAttribute(child,"code",code);
        if (child->FirstChild() && code > 0) {
          utl::getAttribute(child,"type",type,true);
          std::cout <<"\tBodyforce code "<< code;
          if (!type.empty()) std::cout <<" ("<< type <<")";
          VecFunc* f = utl::parseVecFunc(child->FirstChild()->Value(),type);
          if (f) this->setVecProperty(code,Property::BODYLOAD,f);
          std::cout << std::endl;
        }
      }

      else
        return this->Dim::parse(elem);
    }

    if (!mVec.empty())
      this->getIntegrand()->setMaterial(mVec.front());

    return true;
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

protected:
  std::vector<Material*> mVec; //!< Material data

private:
  int aCode; //!< Analytical BC code (used by destructor)
};

#endif
