// $Id$
//==============================================================================
//!
//! \file ModalSim.h
//!
//! \date Aug 29 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for modal analysis of linear dynamics problems.
//!
//==============================================================================

#ifndef _MODAL_SIM_H
#define _MODAL_SIM_H

#include <ios>

class SIMoutput;
class DataExporter;


/*!
  \brief Driver for modal simulation of linear dynamics problems.
  \param[in] infile The input file to parse for time integration setup
  \param[in] nM Number of eigenmodes
  \param[in] dumpModes If \e true, dump projected eigenmode solutions
  \param[in] qstatic If \e true, use quasi-static simulation mode
  \param model The isogeometric finite element model
  \param exporter Result export handler
  \param[in] zero_tol Truncate result values smaller than this to zero
  \param[in] outPrec Number of digits after the decimal point in result print
  \return Exit status
*/

int modalSim (char* infile, size_t nM, bool dumpModes, bool qstatic,
              SIMoutput* model, DataExporter* exporter = nullptr,
              double zero_tol = -1.0, std::streamsize outPrec = 6);

#endif
