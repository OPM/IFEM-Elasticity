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

  //! \brief Parses a data section from an input file.
  //! \details This function allows for specialization of the template
  //! while still reusing as much code as possible.
  //! Only put dimension-specific code in here.
  virtual bool parseDimSpecific(char* keyWord, std::istream& is);

  //! \brief Parses a data section from an XML element.
  //! \details This function allows for specialization of the template
  //! while still reusing as much code as possible.
  //! Only put dimension-specific code in here.
  virtual bool parseDimSpecific(const TiXmlElement* elem);

public:
  //! \brief Returns whether a dual solution is available or not.
  virtual bool haveDualSol() const { return dualS && Dim::dualField; }

private:
  bool dualS; //!< If \e true, also solve the dual problem

  Vector myReact; //!< Nodal reaction forces
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
