// $Id$
//==============================================================================
//!
//! \file NonlinearDriver.h
//!
//! \date Jul 15 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Nonlinear driver for isogeometric finite deformation FEM analysis.
//!
//==============================================================================

#ifndef _NONLINEAR_DRIVER_H
#define _NONLINEAR_DRIVER_H

#include "NonLinSIM.h"
#include "TimeStep.h"

class DataExporter;
class HDF5Restart;
class AdaptiveSetup;


/*!
  \brief Nonlinear driver for isogeometric finite deformation FEM analysis.
  \details This class is derived from NonLinSIM of the IFEM kernel.
  It reimplements the solutionNorms() method to also compute the energy norm
  and other norms of the stress field. In addition, it has the method
  solveProblem() to manage the pseudo-time step loop.
  The driver is also equipped with methods for adaptive mesh refinement based
  on error estimates, through the AdaptiveSetup member.
*/

class NonlinearDriver : public NonLinSIM
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  //! \param sim Reference to the spline FE model
  //! \param[in] linear If \e true, use a linear driver (no Newton iterations)
  //! \param[in] adaptive If \e true, use adaptive mesh refinement
  NonlinearDriver(SIMbase& sim, bool linear, bool adaptive = false);
  //! \brief The destructor deletes the heap-allocated adaptive setup member.
  virtual ~NonlinearDriver();

protected:
  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Computes and prints some solution norm quantities.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] zero_tol Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual bool solutionNorms(const TimeDomain& time,
                             double zero_tol, std::streamsize outPrec);

  //! \brief Prints some solution norm quantities.
  //! \param[in] norm Norm values to print out
  //! \param[in] os The output stream to write the norms to
  virtual void printNorms(const Vector& norm, utl::LogStream& os) const;

  //! \brief Adapts the mesh and restarts solution on new mesh.
  bool adaptMesh(int& aStep);

  //! \brief Calculates and prints out interface force resultants.
  bool calcInterfaceForces(double t);

public:
  //! \brief Reads model data from the specified input file \a *fileName.
  virtual bool read(const char* fileName);

  //! \brief Invokes the main pseudo-time stepping simulation loop.
  //! \param writer HDF5 results exporter
  //! \param restart HDF5 restart handler
  //! \param[in] zero_tol Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  int solveProblem(DataExporter* writer = nullptr,
                   HDF5Restart* restart = nullptr,
                   double zero_tol = 1.0e-8, std::streamsize outPrec = 0)
  {
    return this->solveProblem(writer,restart,nullptr,0,0.0,zero_tol,outPrec);
  }

  //! \brief Invokes the main pseudo-time stepping simulation loop.
  //! \param writer HDF5 results exporter
  //! \param restart HDF5 restart handler
  //! \param oss Output stream for additional ASCII result output
  //! \param[in] printMax Flag for output of maximal stress values
  //! \param[in] dtDump Time increment for dump of ASCII results
  //! \param[in] zero_tol Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  int solveProblem(DataExporter* writer, HDF5Restart* restart,
                   utl::LogStream* oss, bool printMax, double dtDump,
                   double zero_tol = 1.0e-8, std::streamsize outPrec = 0);

  //! \brief Serialize solution state for restarting purposes.
  //! \param data Container for serialized data
  virtual bool serialize(SerializeMap& data) const;
  //! \brief Set solution from a serialized state.
  //! \param[in] data Container for serialized data
  virtual bool deSerialize(const SerializeMap& data);

  //! \brief Accesses the projected solution.
  const Vector& getProjection() const { return proSol.front(); }

  //! \brief Overrides the stop time that was read from the input file.
  void setStopTime(double t) { params.stopTime = t; }

private:
  //! \brief Mark as non-copyable.
  NonlinearDriver(const NonlinearDriver&) = delete;

  //! \brief Mark not comparable.
  bool operator==(const NonlinearDriver&) = delete;

  //! \brief Mark not assignable.
  NonlinearDriver& operator=(const NonlinearDriver&) = delete;

  TimeStep params; //!< Time stepping parameters
  Vectors  proSol; //!< Projected secondary solution
  Matrix   eNorm;  //!< Element norms
  Vectors  gNorm;  //!< Global norms
  char     calcEn; //!< Flag for calculation of solution energy norm
  int      aStep;  //!< Adaptive mesh refinement step
  bool     save0;  //!< If \e true, also save initial configuration

  Vector    myForces;  //!< Interface nodal forces
  Vector    myReacts;  //!< Reaction force container
  RealArray myWeights; //!< Nodal weights for iterface forces

  AdaptiveSetup* adap; //!< Data and methods for adaptive simulation
  std::string inpfile; //!< Model input file, used when adapting mesh
};

#endif
