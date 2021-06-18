// $Id$
//==============================================================================
//!
//! \file SIMLinElModal.h
//!
//! \date Aug 29 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based modal linear elastic FEM analysis.
//!
//==============================================================================

#ifndef _SIM_LIN_EL_MODAL_H
#define _SIM_LIN_EL_MODAL_H

#include "SIMLinEl.h"
#include "SIMmodal.h"


/*!
  \brief Modal driver for isogeometric FEM analysis of linear elastic problems.
  \details The main feature of this class is that it overrides the
  SIMbase::assembleSystem() method, to deal with both assembly of the eigenvalue
  problem, and then also the modal equation system during the time integration.
*/

template<class Dim> class SIMLinElModal : public SIMLinEl<Dim>, public SIMmodal
{
public:
  //! \brief The constructor forwards to the parent class constructors.
  //! \param[in] modes Array of eigenmodes for the elasticity problem
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  SIMLinElModal(std::vector<Mode>& modes, bool checkRHS = false)
    : SIMLinEl<Dim>(nullptr,checkRHS,false), SIMmodal(modes)
  {
    parsed = false;
    alpha1 = alpha2 = 0.0;
  }
  //! \brief Empty destructor.
  virtual ~SIMLinElModal() {}

  using SIMLinEl<Dim>::assembleSystem;
  //! \brief Administers assembly of the linear equation system.
  //! \param[in] time Parameters for time-dependent simulations
  //! \param[in] mSol Previous modal solution vectors
  //!
  //! \details This method assembles the eigenvalue problem in the first call.
  //! Then it is used to build up the modal dynamic equation system during the
  //! time integration. The modal equation system is kept in a separate
  //! AlgEqSystem object (and an associated SAM object), and the pointers to
  //! those modal objects are swapped with the original ones depending on which
  //! system we are currently working on.
  virtual bool assembleSystem(const TimeDomain& time,
                              const Vectors& mSol, bool, bool)
  {
    // Assemble the eigenvalue system
    if (Dim::myProblem->getMode() == SIM::VIBRATION)
      return this->Dim::assembleSystem(TimeDomain(),Vectors());

    if (time.it > 0)
      // Swap back to the full equation system for assembly of load vector
      this->swapSystem(Dim::myEqSys,Dim::mySam);

    else
    {
      // Assemble the load vector of this time step.
      // We need to do this in the first iteration only, as for linear systems
      // the load vector is not supposed to change during the iterations.
      if (!this->Dim::assembleSystem(time,sol,false))
	return false;

      // Extract the load vector in DOF-order
      if (!this->Dim::extractLoadVec(Rhs))
	return false;
    }

    // Assemble the modal equation system
    if (!this->assembleModalSystem(time,mSol,Rhs,
                                   Dim::myProblem->getIntegrationPrm(2),
                                   Dim::myProblem->getIntegrationPrm(3),
                                   alpha1,alpha2))
      return false;

    // Swap the equation systems such that the dynamic simulation driver
    // operates on the modal system
    return this->swapSystem(Dim::myEqSys,Dim::mySam);
  }

  using SIMmodal::expandSolution;
  //! \brief Expands and returns the current dynamic solution.
  //! \param[in] mSol Current modal solution
  //! \param[in] swapBack If \e true, the equation systems are swapped
  virtual const Vectors& expandSolution(const Vectors& mSol, bool swapBack)
  {
    if (!this->expandSolution(mSol,sol))
      sol.clear();

    // Swap back to the full equation system data for postprocessing
    // and assembly of load vector for the next time step
    if (swapBack)
      this->swapSystem(Dim::myEqSys,Dim::mySam);

    return sol;
  }

  //! \brief Returns the current expanded dynamic solution.
  //! \param[in] idx Solution vector index
  virtual const Vector& expandedSolution(int idx) const
  {
    if (idx >= 0 && idx < (int)sol.size())
      return sol[idx];

    static Vector empty;
    return empty;
  }

  //! \brief Returns the number of expanded dynamic solution vectors.
  virtual size_t numExpSolution() const { return sol.size(); }

  //! \brief Serialization support, for the eigenmodes.
  virtual bool serialize(std::map<std::string,std::string>& data) const
  {
    return this->saveModes(data);
  }

  //! \brief Deserialization support (for simulation restart).
  virtual bool deSerialize(const std::map<std::string,std::string>& data)
  {
    return this->restoreModes(data);
  }

  //! \brief Projects the secondary solution associated with the eigenmodes.
  //! \param[out] sesol Control point values of the secondary eigen solutions
  //! \param[out] names Secondary solution component names
  //! \param[in] pMethod Projection method to use
  virtual bool projectModes(Matrices& sesol,
                            std::vector<std::string>& names,
                            SIMoptions::ProjectionMethod pMethod)
  {
    sesol.resize(myModes.size());
    names.resize(Dim::myProblem->getNoFields(2));
    for (size_t c = 0; c < names.size(); c++)
      names[c] = Dim::myProblem->getField2Name(c);

    bool ok = this->setMode(SIM::RECOVERY);
    for (size_t i = 0; i < myModes.size() && ok; i++)
      ok = this->project(sesol[i],myModes[i].eigVec,pMethod);

    return ok;
  }

protected:
  //! \brief Returns the actual integrand.
  //! \details Same as the parent class method, but sets the \a isModal flag.
  virtual Elasticity* getIntegrand()
  {
    if (!Dim::myProblem)
      Dim::myProblem = new LinearElasticity(Dim::dimension,
                                            Elastic::axiSymmetry,
                                            Elastic::GIpointsVTF,true);

    return dynamic_cast<Elasticity*>(Dim::myProblem);
  }

  using SIMLinEl<Dim>::parse;
  //! \brief Parses a data section from an XML element.
  //! \details Overrides the parent class method to do nothing when invoked
  //! during the second time parsing for the time integration setup only.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (parsed)
      IFEM::cout <<"\t(skipped)"<< std::endl;
    else if (!strcasecmp(elem->Value(),"newmarksolver"))
    {
      utl::getAttribute(elem,"alpha1",alpha1);
      utl::getAttribute(elem,"alpha2",alpha2);
    }
    else
      return this->SIMLinEl<Dim>::parse(elem);

    return true;
  }

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details In addition to invoking the inherited method,
  //! this method sets the \a parsed flag, such that the model parsing
  //! is skipped when the input file is parsed for the second time,
  //! for the time integration setup.
  virtual bool preprocessB()
  {
    parsed = true;
    this->setIntegrationPrm(0,alpha1);
    this->setIntegrationPrm(1,alpha2);
    return this->SIMLinEl<Dim>::preprocessB();
  }

private:
  bool   parsed; //!< Set to \e true after the model has been initialized
  double alpha1; //!< Mass-proportional damping parameter
  double alpha2; //!< Stiffness-proportional damping parameter

  Vector  Rhs; //!< Current right-hand-side load vector of the dynamic system
  Vectors sol; //!< Expanded solution vectors from the modal solution
};

#endif
