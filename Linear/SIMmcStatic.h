// $Id$
//==============================================================================
//!
//! \file SIMmcStatic.h
//!
//! \date Feb 12 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Monolithic coupling of linear static simulators.
//!
//==============================================================================

#ifndef _SIM_MC_STATIC_H_
#define _SIM_MC_STATIC_H_

#include "SIMmultiCpl.h"

class DataExporter;


/*!
  \brief Class for monolithic coupled linear static simulators.
*/

class SIMmcStatic : public SIMmultiCpl
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit SIMmcStatic(const std::vector<SIMoutput*>& sims) : SIMmultiCpl(sims)
  {
    myHeading = "Coupled linear static solver";
  }
  //! \brief Empty destructor.
  virtual ~SIMmcStatic() {}

  //! \brief Solves a linear static fully coupled multi-simulator problem.
  int solveStatic(const char* inpfile, DataExporter* exporter,
                  double zero_tol, int outPrec);
};

#endif
