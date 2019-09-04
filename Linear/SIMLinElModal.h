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
    : SIMLinEl<Dim>(checkRHS,false), SIMmodal(modes) {}
  //! \brief Empty destructor.
  virtual ~SIMLinElModal() {}

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

    // Calculate the dynamic solution from the previous modal solution
    if (!this->expandSolution(mSol,sol))
      return false;

    // Swap back the to the full equation system for assembly of load vector.
    // Does nothing in the first call when the modal system is not allocated.
    this->swapSystem(Dim::myEqSys,Dim::mySam);

    if (time.it == 0)
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
                                   Dim::myProblem->getIntegrationPrm(3)))
      return false;

    // Swap the equation systems such that the dynamic simulation driver
    // operates on the modal system
    return this->swapSystem(Dim::myEqSys,Dim::mySam);
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
    if (Dim::myEqSys)
    {
      IFEM::cout <<"\t(skipped)"<< std::endl;
      return true;
    }

    return this->SIMLinEl<Dim>::parse(elem);
  }

private:
  Vector  Rhs; //!< Current right-hand-side load vector of the dynamic system
  Vectors sol; //!< Expanded solution vectors from the modal solution
};

#endif
