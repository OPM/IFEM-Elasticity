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

#include "SIMElasticity.h"
#include "LinearElasticity.h"
#include "ReactionsOnly.h"
#include "AlgEqSystem.h"
#include "SIM2D.h"
#include "SIM3D.h"


/*!
  \brief Driver class for isogeometric FEM analysis of linear elastic problems.
*/

template<class Dim> class SIMLinEl : public SIMElasticity<Dim>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  //! \param[in] ds If \e true, also solve the dual problem
  SIMLinEl(bool checkRHS, bool ds) : SIMElasticity<Dim>(checkRHS), dualS(ds) {}
  //! \brief Constructor for coupled multi-dimensional simulators.
  //! \param[in] head Header identifying this sub-simulator.
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  SIMLinEl(const char* head, bool checkRHS) : SIMElasticity<Dim>(checkRHS)
  {
    Dim::myHeading = head;
    dualS = false;
  }

  //! \brief Empty destructor.
  virtual ~SIMLinEl() {}

  using SIMElasticity<Dim>::solveSystem;
  //! \brief Solves the assembled linear system of equations for a given load.
  //! \param[out] solution Global primary solution vector
  //! \param[in] printSol Print solution if its size is less than \a printSol
  //! \param[out] rCond Reciprocal condition number
  //! \param[in] compName Solution name to be used in norm output
  //! \param[in] newLHS If \e false, reuse the LHS-matrix from previous call
  //! \param[in] idxRHS Index to the right-hand-side vector to solve for
  //!
  //! This overloaded version also computes the reaction forces along a given
  //! boundary. This requires an additional assembly loop calculating the
  //! internal forces only, since we only are doing a linear solve here.
  virtual bool solveSystem(Vector& solution, int printSol, double* rCond,
                           const char* compName, bool newLHS, size_t idxRHS)
  {
    if (!this->Dim::solveSystem(solution,printSol,rCond,compName,newLHS,idxRHS))
      return false;
    else if (idxRHS > 0 || !this->haveBoundaryReactions())
      return true;

    LinearElasticity* prob = dynamic_cast<LinearElasticity*>(Dim::myProblem);
    if (!prob) return true;

    // Assemble the reaction forces. Strictly, we only need to assemble those
    // elements that have nodes on the Dirichlet boundaries, but...
    prob->setReactionIntegral(new ReactionsOnly(myReact,Dim::mySam,Dim::adm));
    AlgEqSystem* tmpEqSys = Dim::myEqSys;
    Dim::myEqSys = nullptr;
    int oldlevel = Dim::msgLevel;
    Dim::msgLevel = 1;
    bool ok = this->setMode(SIM::RHS_ONLY) && this->assembleSystem({solution});

    // Print out the reaction forces
    Vector Rforce;
    if (ok && this->getBoundaryReactions(Rforce))
    {
      IFEM::cout <<"Reaction force     :";
      for (double f : Rforce) IFEM::cout <<" "<< utl::trunc(f);
      IFEM::cout << std::endl;
    }

    Dim::msgLevel = oldlevel;
    Dim::myEqSys = tmpEqSys;
    prob->setReactionIntegral(nullptr);

    return ok;
  }

  //! \brief Returns current reaction force vector.
  virtual const Vector* getReactionForces() const
  {
    return myReact.empty() ? nullptr : &myReact;
  }

  //! \brief Performs static condensation of the linear equation system.
  //! \param[out] Kred Reduced System matrix
  //! \param[out] Rred Associated reduced right-hand-side vector
  virtual bool staticCondensation(Matrix& Kred, Vector& Rred)
  {
    // Assemble [K] and {R}
    this->setMode(SIM::STATIC);
    this->setQuadratureRule(Dim::opt.nGauss[0],true,true);
    if (!this->initSystem(LinAlg::SPARSE))
      return false;
    if (!this->assembleSystem())
      return false;

    return Dim::myEqSys->staticCondensation(Kred,Rred,myRetainNodes,
                                            0,recFile.c_str());
  }

  //! \brief Recovers internal displacements from supernode displacements.
  //! \param[in] Sdisp Local displacements at the supernodes
  //! \param[out] fullDisp Local displacement vector for the entire FE model
  virtual bool recoverInternals(const Vector& Sdisp, Vector& fullDisp)
  {
    if (Rmat.empty() && !AlgEqSystem::readRecoveryMatrix(Rmat,recFile.c_str()))
      return false;

    return Dim::myEqSys->recoverInternals(Rmat,myRetainNodes,Sdisp,fullDisp);
  }

  //! \brief Returns the superelement file name, if any.
  virtual std::string getSupelName() const { return supelName; }

  //! \brief Dumps the supernode coordinates to file.
  virtual void dumpSupernodes(std::ostream& os) const
  {
    os << myRetainNodes.size() << std::endl;
    for (int node : myRetainNodes)
      os << static_cast<Vec3>(this->getNodeCoord(node)) << std::endl;
  }

protected:
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    if (!Dim::myProblem)
      Dim::myProblem = new LinearElasticity(Dim::dimension,
                                            Elastic::axiSymmetry,
                                            Elastic::GIpointsVTF);

    return dynamic_cast<Elasticity*>(Dim::myProblem);
  }

  //! \brief Parses the analytical solution from an input stream.
  virtual bool parseAnaSol(char* keyWord, std::istream& is)
  {
    int code = -1;
    char* cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"EXPRESSION",10))
    {
      IFEM::cout <<"\nAnalytical solution: Expression"<< std::endl;
      int lines = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
      code = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
      if (!Dim::mySol)
        Dim::mySol = new AnaSol(is,lines,false);
    }
    else if (!this->parseDimSpecific(cline))
    {
      std::cerr <<"  ** SIMLinEl::parse: Invalid analytical solution "
                << cline <<" (ignored)"<< std::endl;
      return true;
    }

    // Define the analytical boundary traction field
    if (code == -1)
      code = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
    if (code > 0 && Dim::mySol->getStressSol())
    {
      IFEM::cout <<"Pressure code "<< code
                 <<": Analytical traction"<< std::endl;
      this->setPropertyType(code,Property::NEUMANN);
      Dim::myTracs[code] = new TractionField(*Dim::mySol->getStressSol());
    }

    return true;
  }

  //! \brief Parses the analytical solution from an XML element.
  virtual bool parseAnaSol(const TiXmlElement* elem)
  {
    IFEM::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;

    std::string type;
    utl::getAttribute(elem,"type",type,true);
    if (type == "expression" || type == "fields")
    {
      type[0] = toupper(type[0]);
      IFEM::cout <<"\tAnalytical solution: "<< type << std::endl;
      if (!Dim::mySol)
        Dim::mySol = new AnaSol(elem,false);
    }
    else if (!this->parseDimSpecific(elem,type))
    {
      std::cerr <<"  ** SIMLinEl::parse: Invalid analytical solution "
                << type <<" (ignored)"<< std::endl;
      return true;
    }

    // Define the analytical boundary traction field
    int code = 0;
    utl::getAttribute(elem,"code",code);
    if (code > 0 && Dim::mySol && Dim::mySol->getStressSol())
    {
      IFEM::cout <<"\tNeumann code "<< code
                 <<": Analytical traction"<< std::endl;
      this->setPropertyType(code,Property::NEUMANN);
      Dim::myTracs[code] = new TractionField(*Dim::mySol->getStressSol());
    }

    return true;
  }

  using SIMElasticity<Dim>::parse;
  //! \brief Parses a data section from an XML element
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem)
  {
    // Lambda function for parsing the static condensation tag.
    auto&& parseSC = [this](const TiXmlElement* elem)
    {
      IFEM::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;
      if (utl::getAttribute(elem,"supelName",supelName))
      {
        recFile = supelName;
        size_t idot = recFile.find_last_of('.');
        if (idot < recFile.size())
          recFile.replace(idot,std::string::npos,".rec");
        else
          recFile.append(".rec");
      }

      std::string rset;
      const TiXmlElement* child = elem->FirstChildElement();
      for (; child; child = child->NextSiblingElement())
        if (!strcasecmp(child->Value(),"retained"))
          if (utl::getAttribute(child,"set",rset))
          {
            const TopEntity& retEnt = this->getEntity(rset);
            if (retEnt.empty())
              IFEM::cout <<"  ** Undefined topology set "<< rset << std::endl;
            else
            {
              myRetainSet.insert(retEnt.begin(),retEnt.end());
              IFEM::cout <<"\tRetained topology set: "<< rset << std::endl;
            }
          }
      return true;
    };

    if (!this->SIMElasticity<Dim>::parse(elem))
      return false;
    else if (!strcasecmp(elem->Value(),"staticCondensation"))
      return parseSC(elem);
    else if (strcasecmp(elem->Value(),SIMElasticity<Dim>::myContext.c_str()))
      return true;

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"staticCondensation"))
        return parseSC(child);

    return true;
  }

  //! \brief Parses the analytical solution from an input stream.
  //! \details This function allows for specialization of the template
  //! while still reusing as much code as possible.
  //! Only put dimension-specific code in here.
  bool parseDimSpecific(char* cline);

  //! \brief Parses the analytical solution from an XML element.
  //! \details This function allows for specialization of the template
  //! while still reusing as much code as possible.
  //! Only put dimension-specific code in here.
  bool parseDimSpecific(const TiXmlElement* elem, const std::string& type);

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to resolve the topology sets
  //! into node numbers specified for static condensation.
  virtual bool preprocessB()
  {
    bool ok = SIMElasticity<Dim>::preprocessB();
    if (!ok) return false;

    for (const TopItem& titem : myRetainSet)
      if (!this->getTopItemNodes(titem,myRetainNodes))
      {
        IFEM::cout <<" *** SIMLinEl::preprocessB: Invalid topology set "
                   << titem << std::endl;
        ok = false;
      }

    if (!ok || myRetainNodes.empty())
      return ok;

    std::sort(myRetainNodes.begin(),myRetainNodes.end());
    IFEM::cout <<"\nRetained nodes in static condensation:";
    for (size_t i = 0; i < myRetainNodes.size(); i++)
      IFEM::cout << (i%10 ? " " : "\n") << myRetainNodes[i];
    IFEM::cout <<"\n"<< std::endl;
    return true;
  }

public:
  //! \brief Returns whether a dual solution is available or not.
  virtual bool haveDualSol() const { return dualS && Dim::dualField; }

private:
  bool dualS; //!< If \e true, also solve the dual problem

  Vector myReact; //!< Nodal reaction forces

  TopEntity myRetainSet;   //!< Topology set for the retained DOFs
  IntVec    myRetainNodes; //!< List of retained nodes in static condensation

  std::string supelName; //!< Name of superelement file
  std::string recFile;   //!< Name of displacement recovery file

  Matrix Rmat; //!< Static condensation recovery matrix
};

typedef SIMLinEl<SIM2D> SIMLinEl2D; //!< 2D specific driver
typedef SIMLinEl<SIM3D> SIMLinEl3D; //!< 3D specific driver

//! \brief Template specialization - 2D specific input parsing.
template<> bool SIMLinEl2D::parseDimSpecific(char* cline);
//! \brief Template specialization - 2D specific input parsing.
template<> bool SIMLinEl2D::parseDimSpecific(const TiXmlElement* elem,
                                             const std::string& type);

//! \brief Template specialization - 3D specific input parsing.
template<> bool SIMLinEl3D::parseDimSpecific(char* cline);
//! \brief Template specialization - 3D specific input parsing.
template<> bool SIMLinEl3D::parseDimSpecific(const TiXmlElement* elem,
                                             const std::string& type);

#endif
