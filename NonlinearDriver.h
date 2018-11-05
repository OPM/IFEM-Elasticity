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


/*!
  \brief Nonlinear driver for isogeometric finite deformation FEM analysis.
  \details This class is derived from NonLinSIM of the IFEM kernel.
  It reimplements the \a solutionNorms method to also compute the energy norm
  and other norms of the stress field. In addition, it has the method
  \a solveProblem to manage the pseudo-time step loop.
*/

class NonlinearDriver : public NonLinSIM
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  //! \param sim Reference to the spline FE model
  //! \param linear If \e true, use a linear driver (no Newton iterations)
  explicit NonlinearDriver(SIMbase& sim, bool linear = false);
  //! \brief Empty destructor.
  virtual ~NonlinearDriver() {}

protected:
  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

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

public:
  //! \brief Invokes the main pseudo-time stepping simulation loop.
  //! \param writer HDF5 results exporter
  //! \param restart HDF5 restart handler
  //! \param oss Output stream for additional ASCII result output
  //! \param[in] dtDump Time increment for dump of ASCII results
  //! \param[in] zero_tol Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  int solveProblem(DataExporter* writer, HDF5Restart* restart, utl::LogStream* oss,
                   double dtDump, double zero_tol, std::streamsize outPrec);

  //! \brief Serialize solution state for restarting purposes.
  //! \param data Container for serialized data
  virtual bool serialize(SerializeMap& data) const;
  //! \brief Set solution from a serialized state.
  //! \param[in] data Container for serialized data
  virtual bool deSerialize(const SerializeMap& data);

  //! \brief Accesses the projected solution.
  const Vector& getProjection() const { return proSol; }

  //! \brief Overrides the stop time that was read from the input file.
  void setStopTime(double t) { params.stopTime = t; }
  //! \brief Flag calculation of solution energy norms.
  void calculateEnergy(char flag) { calcEn = flag; }
  //! \brief Flag that we are doing a linear analysis only.
  void setLinear() { iteNorm = NONE; }

private:
  TimeStep params; //!< Time stepping parameters
  char     calcEn; //!< Flag for calculation of solution energy norm
  Matrix   proSol; //!< Projected secondary solution
};

#endif
