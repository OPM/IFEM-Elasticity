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

  //! \brief Serialize solution state for restarting purposes.
  //! \param data Container for serialized data
  virtual bool serialize(HDF5Restart::SerializeData& data) const
  {
    return (model.serialize(data) &&
            this->NewmarkDriver<NewmarkSIM>::serialize(data));
  }

  //! \brief Set solution from a serialized state.
  //! \param[in] data Container for serialized data
  virtual bool deSerialize(const HDF5Restart::SerializeData& data)
  {
    return (model.deSerialize(data) &&
            this->NewmarkDriver<NewmarkSIM>::deSerialize(data));
  }
};


int modalSim (char* infile, size_t nM, SIMoutput* model, DataExporter* exporter)
{
  ModalDriver simulator(*model);

  // Read time integration setup
  if (!strcasestr(infile,".xinp") || !simulator.readXML(infile))
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
  int status = simulator.solveProblem(exporter,restart);

  delete restart;
  return status;
}
