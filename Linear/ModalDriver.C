// $Id$
//==============================================================================
//!
//! \file ModalDriver.C
//!
//! \date Aug 29 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for modal analysis of linear dynamics problems.
//!
//==============================================================================

#include "NewmarkDriver.h"
#include "NewmarkSIM.h"
#include "SIMmodal.h"


/*!
  \brief Driver for modal analysis of linear dynamic problems.
*/

class ModalDriver : public NewmarkDriver<NewmarkSIM>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  //! \param sim Reference to the spline FE model
  explicit ModalDriver(SIMbase& sim) : NewmarkDriver<NewmarkSIM>(sim) {}
  //! \brief Empty destructor.
  virtual ~ModalDriver() {}

  //! \brief Returns the number of solution vectors.
  virtual size_t numSolution() const
  {
    return dynamic_cast<SIMmodal*>(&model)->numExpSolution();
  }

  //! \brief Calculates the current real solution vectors.
  virtual const Vectors& realSolutions()
  {
    return dynamic_cast<SIMmodal*>(&model)->expandSolution(solution,true);
  }

  //! \brief Returns a const reference to the current real solution vector.
  virtual const Vector& realSolution(int i = 0) const
  {
    return dynamic_cast<SIMmodal*>(&model)->expandedSolution(i);
  }
};


int modalSim (char* infile, size_t nM, SIMoutput* model, DataExporter* exporter)
{
  ModalDriver simulator(*model);

  // Read time integration setup
  if (!simulator.read(infile))
    return 1;

  // Initialize the modal time-domain simulation
  simulator.initPrm();
  simulator.initSolution(nM,3);
  simulator.printProblem();

  // Save FE model to VTF file for visualization
  if (model->opt.format >= 0)
    if (!simulator.saveModel(infile))
      return 2;

  // Initialize the linear equation system.
  // Actually, we don't need any system matrices here
  // since we are only integrating the external load vector in space.
  if (!model->initSystem(LinAlg::DENSE,0,1))
    return 3;

  // Run the modal time integration
  return simulator.solveProblem(exporter,nullptr);
}
