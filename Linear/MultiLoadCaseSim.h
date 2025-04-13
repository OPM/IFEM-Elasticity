// $Id$
//==============================================================================
//!
//! \file MultiLoadCaseSim.h
//!
//! \date May 05 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for linear multi-load-case static problems.
//!
//==============================================================================

#ifndef _MULTI_LOAD_CASE_SIM_H
#define _MULTI_LOAD_CASE_SIM_H

#include <ios>

class SIMoutput;


/*!
  \brief Driver for modal simulation of linear dynamics problems.
  \param[in] infile The input file to parse for load case setup
  \param model The isogeometric finite element model
  \param[in] fixDup If \e true, merge duplicated FE nodes
  \param[in] dumpNodeMap If \e true, export node mapping to HDF5 file
  \param[in] zero_tol Truncate result values smaller than this to zero
  \param[in] outPrec Number of digits after the decimal point in result print
  \return Exit status
*/

int mlcSim (char* infile, SIMoutput* model,
            bool fixDup = false, bool dumpNodeMap = false,
            double zero_tol = 1.0e-8, std::streamsize outPrec = 6);

#endif
