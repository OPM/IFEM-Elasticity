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
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Returns the actual integrand.
  virtual ElasticBase* getIntegrand();

  //! \brief Performs some pre-processing tasks on the FE model.
  virtual bool preprocessB();

  //! \brief Creates a property set for contact condition on an entity set.
  //! \param[in] slaveSet Name of the slave boundary entity set
  //! \param[out] code Property code associated with the contact set
  virtual bool createContactSet(const std::string& slaveSet, int& code);

  //! \brief Specialized preprocessing performed before assembly initialization.
  //! \param ngnod Total number of nodal points in this model
  virtual bool preprocessBeforeAsmInit(int& ngnod);

  //! \brief Assembles problem-dependent discrete terms, if any.
  virtual bool assembleDiscreteTerms(const IntegrandBase* problem,
                                     const TimeDomain&);

public:
  //! \brief Prints out problem-specific data to the log stream.
  virtual bool printProblem() const;

  //! \brief Updates the time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] time Current time
  //! \param[in] prevSol Pointer to previous primary solution in DOF-order
  virtual bool updateDirichlet(double time, const Vector* prevSol);

  //! \brief Updates Problem-dependent state based on current solution.
  //! \param[in] solution Current primary solution vector
  virtual bool updateConfiguration(const Vector& solution);

  //! \brief Prints interface force resultants associated with given boundaries.
  //! \param[in] sf Internal nodal forces
  //! \param weights Nodal weights (in case some nodes are present in more sets)
  virtual void printIFforces(const Vector& sf, RealArray& weights);

  //! \brief Writes current model geometry to the VTF-file.
  //! \param nBlock Running result block counter
  //! \param[in] inpFile File name used to construct the VTF-file name from
  //!
  //! \details This method is overridden to also write out the contact bodies.
  virtual bool writeGlvG(int& nBlock, const char* inpFile, bool = true);
  //! \brief Writes contact body movements to the VTF-file.
  //! \param nBlock Running result block counter
  //! \param[in] iStep Load/time step identifier
  virtual bool writeGlvA(int& nBlock, int iStep, double, int) const;

  //! \brief Dumps additional problem-specific results in ASCII format.
  //! \param os Output stream to write the solution data to
  //! \param[in] prec Number of digits after the decimal point
  virtual void dumpMoreResults(double, utl::LogStream& os,
                               std::streamsize prec) const;

private:
  NLoptions nlo; //!< Input options defining the nonlinear formulation
};

#endif
