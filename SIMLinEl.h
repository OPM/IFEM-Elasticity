// $Id$
//==============================================================================
//!
//! \file SIMLinEl.h
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#ifndef _SIM_LIN_EL_H
#define _SIM_LIN_EL_H

#include "LinearElasticity.h"
#include "LinIsotropic.h"
#include "AnaSol.h"
#include "Functions.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "Utilities.h"
#include "tinyxml.h"


/*!
  \brief Driver class for isogeometric FEM analysis of elasticity problems.
  \details The class incapsulates data and methods for solving linear elasticity
  problems using NURBS-based finite elements. It reimplements the parse method
  of the parent class.
*/

template<class Dim> class SIMLinEl : public Dim
{
public:
  //! \brief Default constructor.
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  SIMLinEl(bool checkRHS = false) : Dim(Dim::dimension,0,checkRHS) { aCode = 0; }

  //! \brief The destructor frees the dynamically allocated material properties.
  virtual ~SIMLinEl()
  {
    // To prevent the SIMbase destructor try to delete already deleted functions
    if (aCode > 0)
      Dim::myVectors.erase(aCode);

    for (size_t i = 0; i < mVec.size(); i++)
      delete mVec[i];
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

    this->SIMbase::clearProperties();
  }

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented inserting a call to \a getIntegrand.
  //! This makes sure the integrand has been allocated in case of minimum input.
  //! It also resolves inhomogeneous boundary condition fields in case they are
  //! derived from the analytical solution.
  virtual bool preprocess(const std::vector<int>& ignored, bool fixDup)
  {
    if (!Dim::myProblem)
    {
      this->getIntegrand();
      if (Dim::myPid == 0)
        this->printProblem(std::cout);
    }

    if (Dim::mySol) // Define analytical boundary condition fields
      for (PropertyVec::iterator p = Dim::myProps.begin();
                                 p != Dim::myProps.end(); p++)
        if (p->pcode == Property::DIRICHLET_ANASOL)
        {
          if (!Dim::mySol->getVectorSol())
            p->pcode = Property::UNDEFINED;
          else if (aCode == abs(p->pindx))
            p->pcode = Property::DIRICHLET_INHOM;
          else if (aCode == 0)
          {
            aCode = abs(p->pindx);
            Dim::myVectors[aCode] = Dim::mySol->getVectorSol();
            p->pcode = Property::DIRICHLET_INHOM;
          }
          else
            p->pcode = Property::UNDEFINED;
        }
        else if (p->pcode == Property::NEUMANN_ANASOL)
        {
          if (Dim::mySol->getStressSol())
          {
            p->pcode = Property::NEUMANN;
            Dim::myTracs[p->pindx] = new TractionField(*Dim::mySol->getStressSol());
          }
          else
            p->pcode = Property::UNDEFINED;
        }

    return this->Dim::preprocess(ignored,fixDup);
  }

  static bool planeStrain; //!< Plane strain/stress option - 2D only
  static bool axiSymmetry; //!< Axisymmtry option - 2D only

private:
  //! \brief Returns the actual integrand.
  Elasticity* getIntegrand()
  {
    if (Dim::myProblem)
      return dynamic_cast<Elasticity*>(Dim::myProblem);

    Elasticity* elp = new LinearElasticity(Dim::dimension,
                                         Dim::dimension == 2?axiSymmetry:false);
    Dim::myProblem = elp;

    return elp;
  }

protected:
  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is)
  {
    char* cline = 0;
    int nmat = 0;
    int nConstPress = 0;
    int nLinearPress = 0;

    if (parseDimSpecific(keyWord, is))
      return true;
    else if (!strncasecmp(keyWord,"ISOTROPIC",9))
    {
      nmat = atoi(keyWord+10);
      if (Dim::myPid == 0)
        std::cout <<"\nNumber of isotropic materials: "<< nmat << std::endl;

      for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
      {
        int code = atoi(strtok(cline," "));
        if (code > 0)
          this->setPropertyType(code,Property::MATERIAL,mVec.size());

        double E   = atof(strtok(NULL," "));
        double nu  = atof(strtok(NULL," "));
        double rho = atof(strtok(NULL," "));
        if (Dim::dimension == 3)
          mVec.push_back(new LinIsotropic(E,nu,rho));
        if (Dim::dimension == 2)
          mVec.push_back(new LinIsotropic(E,nu,rho,!planeStrain,axiSymmetry));
        if (Dim::myPid == 0)
          std::cout <<"\tMaterial code "<< code <<": "
                    << E <<" "<< nu <<" "<< rho << std::endl;
      }
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
      press.ldim = Dim::dimension-1;

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
        if (press.lindx < 1 || press.lindx > 4+2*(Dim::dimension-2))
        {
          std::cerr <<" *** SIMLinEl3D::parse: Invalid face index "
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
            Dim::myTracs[1+i] = new PressureField(utl::parseRealFunc(cline,p),pdir);
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
            if (Dim::dimension == 3)
              mVec.push_back(new LinIsotropic(E,nu,rho));
            if (Dim::dimension == 2)
              mVec.push_back(new LinIsotropic(E,nu,rho,
                                              !planeStrain,axiSymmetry));
          }
          else
          {
            int patch = atoi(cline);
            int pid = this->getLocalPatchIndex(patch);
            if (pid < 0) return false;
            if (pid < 1) continue;

            std::cout <<"\tMaterial for P"<< patch
                      <<": "<< E <<" "<< nu <<" "<< rho << std::endl;
            Dim::myProps.push_back(Property(Property::MATERIAL,mVec.size(),pid,3));
            if (Dim::dimension == 3)
              mVec.push_back(new LinIsotropic(E,nu,rho));
            if (Dim::dimension == 2)
              mVec.push_back(new LinIsotropic(E,nu,rho,
                                              !planeStrain,axiSymmetry));
          }
      }
    }

    else
      return this->Dim::parse(keyWord,is);

    int npres = nConstPress + nLinearPress;
    if (npres > 0)
    {
      if (Dim::myPid == 0)
        std::cout <<"\nNumber of pressures: "<< npres << std::endl;
      for (int i = 0; i < npres && (cline = utl::readLine(is)); i++)
      {
        int code = atoi(strtok(cline," "));
        int pdir = atoi(strtok(NULL," "));
        double p = atof(strtok(NULL," "));
        if (Dim::myPid == 0)
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
      if (parseDimSpecific(child))
        continue;

      if (!strcasecmp(child->Value(),"isotropic")) {
        int code = this->parseMaterialSet(child,mVec.size());

        double E = 1000.0, nu = 0.3, rho = 1.0;
        utl::getAttribute(child,"E",E);
        utl::getAttribute(child,"nu",nu);
        utl::getAttribute(child,"rho",rho);

        if (Dim::dimension == 3)
          mVec.push_back(new LinIsotropic(E,nu,rho));
        if (Dim::dimension == 2)
          mVec.push_back(new LinIsotropic(E,nu,rho, !planeStrain,axiSymmetry));
        if (Dim::myPid == 0)
          std::cout <<"\tMaterial code "<< code <<": "
                    << E <<" "<< nu <<" "<< rho << std::endl;
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
      } else
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

  //! \brief Parse a data section from an input file
  //! \details This function allows for specialization of the template
  //!          while still reusing as much code as possible. Only put
  //!          dimension-specific code in here
  bool parseDimSpecific(char* keyWord, std::istream& is);

  //! \brief Parse a data section from an input file
  //! \details This function allows for specialization of the template
  //!          while still reusing as much code as possible. Only put
  //!          dimension-specific code in here
  bool parseDimSpecific(const TiXmlElement* elem);

  //! \brief Initializes the body load properties for current patch.
  //! \param[in] patchInd 1-based patch index
  virtual bool initBodyLoad(size_t patchInd)
  {
    Elasticity* elp = dynamic_cast<Elasticity*>(Dim::myProblem);
    if (!elp) return false;

    typename Dim::VecFuncMap::const_iterator it = Dim::myVectors.end();
    for (size_t i = 0; i < Dim::myProps.size(); i++)
      if (Dim::myProps[i].pcode == Property::BODYLOAD &&
          Dim::myProps[i].patch == patchInd)
        if ((it = Dim::myVectors.find(Dim::myProps[i].pindx)) != Dim::myVectors.end())
          break;

    elp->setBodyForce(it == Dim::myVectors.end() ? NULL : it->second);
    return true;
  }

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd)
  {
    Elasticity* elp = dynamic_cast<Elasticity*>(Dim::myProblem);
    if (!elp) return false;

    typename Dim::VecFuncMap::const_iterator  vit = Dim::myVectors.find(propInd);
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

typedef SIMLinEl<SIM2D> SIMLinEl2D; //!< 2D specific driver
typedef SIMLinEl<SIM3D> SIMLinEl3D; //!< 3D specific driver

//! \brief Template specialization - 2D specific input parsing.
template<> bool SIMLinEl2D::parseDimSpecific(char* keyWord, std::istream& is);
//! \brief Template specialization - 2D specific input parsing.
template<> bool SIMLinEl2D::parseDimSpecific(const TiXmlElement* elem);

//! \brief Template specialization - 3D specific input parsing.
template<> bool SIMLinEl3D::parseDimSpecific(char* keyWord, std::istream& is);
//! \brief Template specialization - 3D specific input parsing.
template<> bool SIMLinEl3D::parseDimSpecific(const TiXmlElement* elem);

#endif
