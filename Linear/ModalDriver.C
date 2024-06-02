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

#include "ModalDriver.h"
#include "SIMmodal.h"

using Parent = NewmarkDriver<NewmarkSIM>; //!< Convenience renaming


const Vectors& ModalDriver::realSolutions ()
{
  return dynamic_cast<SIMmodal*>(&model)->expandSolution(solution,true);
}


const Vector& ModalDriver::realSolution (int i) const
{
  return dynamic_cast<SIMmodal*>(&model)->expandedSolution(i);
}


size_t ModalDriver::numSolution () const
{
  return dynamic_cast<SIMmodal*>(&model)->numExpSolution();
}


bool ModalDriver::serialize (HDF5Restart::SerializeData& data) const
{
  return model.serialize(data) && this->Parent::serialize(data);
}


bool ModalDriver::deSerialize (const HDF5Restart::SerializeData& data)
{
  return model.deSerialize(data) && this->Parent::deSerialize(data);
}


void ModalDriver::dumpResults (double time, utl::LogStream& os,
                               std::streamsize precision, bool formatted) const
{
  model.dumpResults(this->realSolution(),time,os,formatted,precision);
}


void ModalDriver::dumpModes (utl::LogStream& os,
                             std::streamsize precision) const
{
  // Project the secondary eigensolution
  SIMoptions::ProjectionMap::const_iterator pit = model.opt.project.begin();
  if (pit == model.opt.project.end()) return; // No projection specified

  Matrices sesol;
  std::vector<std::string> names;
  if (!dynamic_cast<SIMmodal*>(&model)->projectModes(sesol,names,pit->first))
    return; // Projection failure (ignored).

  // Formatted output, use scientific notation with fixed field width
  std::streamsize flWidth = 8 + precision;
  std::streamsize oldPrec = os.precision(precision);
  std::ios::fmtflags oldF = os.flags(std::ios::scientific | std::ios::right);
  int modeNo = 0;
  for (const Matrix& ssol : sesol)
  {
    os <<"\n\n >>> Secondary solution for mode "<< ++modeNo <<":\n\nCtrl.pt.";
    for (const std::string& cmp : names) os << std::setw(flWidth) << cmp;
    for (size_t ipt = 1; ipt <= ssol.cols(); ipt++)
    {
      os <<"\n"<< std::setw(7) << ipt <<" ";
      for (size_t cmp = 1; cmp <= ssol.rows(); cmp++)
        os << std::setw(flWidth) << utl::trunc(ssol(cmp,ipt));
    }
  }
  os <<"\n"<< std::endl;
  os.precision(oldPrec);
  os.flags(oldF);
}


bool ModalDriver::predictStep (TimeStep& tp)
{
  return qstatic ? true : this->Parent::predictStep(tp);
}


bool ModalDriver::correctStep (TimeStep& tp, bool)
{
  if (qstatic)
  {
    solution.front() = linsol;
    model.updateRotations(solution.front());
    return model.updateConfiguration(solution.front());
  }
  else
    return this->Parent::correctStep(tp);
}


SIM::ConvStatus ModalDriver::checkConvergence (TimeStep& tp)
{
  return qstatic ? SIM::CONVERGED : this->Parent::checkConvergence(tp);
}


int modalSim (char* infile, size_t nM, bool dumpModes, bool qstatic,
              SIMoutput* model, DataExporter* exporter,
              double zero_tol, std::streamsize outPrec)
{
  ModalDriver simulator(*model,qstatic);

  // Print out control point stresses for the eigenmodes
  if (dumpModes)
    simulator.dumpModes(IFEM::cout,10);

  // Read time integration setup
  if (!strcasestr(infile,".xinp") || !simulator.readXML(infile))
    return 1;

  // Initialize the modal time-domain simulation
  simulator.initPrm();
  // Set beta=0 for quasi-static simulation
  if (qstatic) model->setIntegrationPrm(2,0.0);
  simulator.initSolution(nM, qstatic ? 1 : 3);
  simulator.printProblem();

  // Save FE model to VTF file for visualization
  if (model->opt.format >= 0 && !simulator.saveModel(infile))
    return 2;

  // Initialize the linear equation system.
  // Actually, we don't need any system matrices here
  // since we are only integrating the external load vector in space.
  if (!model->initSystem(LinAlg::DENSE,0,1))
    return 3;

  // Load solution state from serialized data in case of restart
  if (!simulator.checkForRestart())
    return 4;

  // Check for restart output
  HDF5Restart* restart = nullptr;
  if (model->opt.restartInc > 0)
  {
    std::string hdf5file(infile);
    if (!model->opt.hdf5.empty())
      hdf5file = model->opt.hdf5 + "_restart";
    else
      hdf5file.replace(hdf5file.find_last_of('.'),std::string::npos,"_restart");
    IFEM::cout <<"\nWriting HDF5 file "<< hdf5file <<".hdf5"<< std::endl;
    restart = new HDF5Restart(hdf5file,model->getProcessAdm(),
                              model->opt.restartInc);
  }

  // Run the modal time integration
  int status = simulator.solveProblem(exporter,restart,zero_tol,outPrec);
  delete restart;
  return status;
}
