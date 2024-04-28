// $Id$
//==============================================================================
//!
//! \file SIMLinKLModal.h
//!
//! \date Nov 21 2023
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for modal linear Kirchhoff-Love FEM analysis.
//!
//==============================================================================

#ifndef _SIM_LIN_KL_MODAL_H
#define _SIM_LIN_KL_MODAL_H

#include "SIMLinElKL.h"
#include "SIMmodal.h"


/*!
  \brief Modal driver for isogeometric FEM analysis of Kirchhoff-Love problems.
  \details The main feature of this class is that it overrides the
  SIMbase::assembleSystem() method, to deal with both assembly of the eigenvalue
  problem, and then also the modal equation system during the time integration.
*/

class SIMLinKLModal : public SIMLinElKL, public SIMmodal
{
public:
  //! \brief The constructor forwards to the parent class constructors.
  //! \param[in] modes Array of eigenmodes for the elasticity problem
  //! \param[in] shell If \e true, 3 DOFs per node is assumed, otherwise 1
  explicit SIMLinKLModal(std::vector<Mode>& modes, bool shell = false);
  //! \brief Empty destructor.
  virtual ~SIMLinKLModal() {}

  using SIMLinElKL::assembleSystem;
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
                              const Vectors& mSol, bool, bool);

  using SIMmodal::expandSolution;
  //! \brief Expands and returns the current dynamic solution.
  //! \param[in] mSol Current modal solution
  //! \param[in] swapBck If \e true, the equation systems are swapped
  virtual const Vectors& expandSolution(const Vectors& mSol, bool swapBck);

  //! \brief Returns the current expanded dynamic solution.
  //! \param[in] idx Solution vector index
  virtual const Vector& expandedSolution(int idx) const;

  //! \brief Returns the number of expanded dynamic solution vectors.
  virtual size_t numExpSolution() const { return sol.size(); }

  //! \brief Serialization support, for the eigenmodes.
  virtual bool serialize(std::map<std::string,std::string>& data) const;
  //! \brief Deserialization support (for simulation restart).
  virtual bool deSerialize(const std::map<std::string,std::string>& data);

  //! \brief Projects the secondary solution associated with the eigenmodes.
  //! \param[out] sesol Control point values of the secondary eigen solutions
  //! \param[out] names Secondary solution component names
  //! \param[in] pMethod Projection method to use
  virtual bool projectModes(Matrices& sesol,
                            std::vector<std::string>& names,
                            SIMoptions::ProjectionMethod pMethod);

protected:
  using SIMLinElKL::parse;
  //! \brief Parses a data section from an XML element.
  //! \details Overrides the parent class method to do nothing when invoked
  //! during the second time parsing for the time integration setup only.
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details In addition to invoking the inherited method,
  //! this method sets the \a parsed flag, such that the model parsing
  //! is skipped when the input file is parsed for the second time,
  //! for the time integration setup.
  virtual bool preprocessB();

private:
  bool   parsed; //!< Set to \e true after the model has been initialized
  double alpha1; //!< Mass-proportional damping parameter
  double alpha2; //!< Stiffness-proportional damping parameter

  Vector  Rhs; //!< Current right-hand-side load vector of the dynamic system
  Vectors sol; //!< Expanded solution vectors from the modal solution
};

#endif
