// $Id$
//==============================================================================
//!
//! \file SIMFiniteDefEl.h
//!
//! \date Dec 18 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution drivers for NURBS-based finite deformation analysis.
//!
//==============================================================================

#ifndef _SIM_FINITE_DEF_EL_H
#define _SIM_FINITE_DEF_EL_H

#include "SIMElasticity.h"
#include "SIMContact.h"
#include "NLoptions.h"


/*!
  \brief Driver class for isogeometric finite deformation analysis.
*/

template<class Dim>
class SIMFiniteDefEl : public SIMElasticity<Dim>, private SIMContact
{
  std::vector<Material*>& mDat; //!< Convenience reference to material data

public:
  //! \brief The constructor initializes the formulation options.
  //! \param[in] rhs If \e true, ensure the model is in a right-hand system
  //! \param[in] options Additional solution formulation options
  SIMFiniteDefEl(bool rhs, const std::vector<int>& options);

protected:
  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  bool parse(char* keyWord, std::istream& is) override;
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  bool parse(const tinyxml2::XMLElement* elem) override;

  //! \brief Parses material data from an XML element.
  //! \param[in] elem The XML element to parse
  virtual Material* parseMaterial(const tinyxml2::XMLElement* elem);

  //! \brief Returns the actual integrand.
  ElasticBase* getIntegrand() override;

  //! \brief Performs some pre-processing tasks on the FE model.
  bool preprocessB() override;

  //! \brief Creates a property set for contact condition on an entity set.
  //! \param[in] slaveSet Name of the slave boundary entity set
  //! \param[out] code Property code associated with the contact set
  bool createContactSet(const std::string& slaveSet, int& code) override;

  //! \brief Specialized preprocessing performed before assembly initialization.
  //! \param ngnod Total number of nodal points in this model
  bool preprocessBeforeAsmInit(int& ngnod) override;

  //! \brief Assembles problem-dependent discrete terms, if any.
  bool assembleDiscreteTerms(const IntegrandBase* problem,
                             const TimeDomain&) override;

public:
  //! \brief Prints out problem-specific data to the log stream.
  bool printProblem() const override;

  //! \brief Updates the time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] time Current time
  //! \param[in] prevSol Pointer to previous primary solution in DOF-order
  //! \param[in] tangent If \e true, use time-derivatives of prescribed values
  bool updateDirichlet(double time, const Vector* prevSol, bool tangent) override;

  //! \brief Updates Problem-dependent state based on current solution.
  //! \param[in] solution Current primary solution vector
  bool updateConfiguration(const Vector& solution) override;

  //! \brief Prints interface force resultants associated with given boundaries.
  //! \param[in] sf Internal nodal forces
  //! \param weights Nodal weights (in case some nodes are present in more sets)
  void printIFforces(const Vector& sf, RealArray& weights) override;

  //! \brief Writes current model geometry to the VTF-file.
  //! \param nBlock Running result block counter
  //! \param[in] time The time from which this (new) geometry applies
  //! \param[in] append If \e true, append new blocks to existing ones, if any
  //!
  //! \details This method is overridden to also write out the contact bodies.
  bool writeGlvG(int& nBlock, double time, bool append) override;
  //! \brief Writes contact body movements to the VTF-file.
  //! \param nBlock Running result block counter
  //! \param[in] iStep Load/time step identifier
  bool writeGlvA(int& nBlock, int iStep, double, int) const override;

  //! \brief Dumps additional problem-specific results in ASCII format.
  //! \param os Output stream to write the solution data to
  //! \param[in] prec Number of digits after the decimal point
  void dumpMoreResults(double, utl::LogStream& os,
                       std::streamsize prec) const override;

private:
  NLoptions nlo; //!< Input options defining the nonlinear formulation
};

#endif
