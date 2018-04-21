// $Id$
//==============================================================================
//!
//! \file SIMLinElBeamC1.h
//!
//! \date Sep 16 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based FE analysis of C1-continuous beams.
//!
//==============================================================================

#ifndef _SIM_LIN_EL_BEAM_C1_H
#define _SIM_LIN_EL_BEAM_C1_H

#include "SIM1D.h"

class Material;


/*!
  \brief Driver class for isogeometric FEM analysis of C1-continuous beams.
*/

class SIMLinElBeamC1 : public SIM1D
{
  /*!
    \brief Struct defining a nodal point load.
  */
  struct PointLoad
  {
    size_t patch; //!< Patch index [0,nPatch>
    int    inod;  //!< Local node number of the closest node
    double xi;    //!< Parameter of the point
    Vec3   X;     //!< Spatial coordinates of the point
    double pload; //!< Load magnitude
    // \brief Default constructor.
    PointLoad() : patch(0), inod(0) { xi = pload = 0.0; }
  };

  typedef std::vector<PointLoad> PloadVec; //!< Point load container

public:
  //! \brief Default constructor.
  SIMLinElBeamC1() : SIM1D(1) {}
  //! \brief Empty destructor.
  virtual ~SIMLinElBeamC1() {}

protected:
  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Performs some pre-processing tasks on the FE model.
  virtual bool preprocessB();

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd);
  //! \brief Initializes the body load properties for current patch.
  //! \param[in] patchInd 1-based patch index
  virtual bool initBodyLoad(size_t patchInd);

  //! \brief Assembles the nodal point loads, if any.
  virtual bool assembleDiscreteTerms(const IntegrandBase*, const TimeDomain&);

  //! \brief Computes problem-dependent external energy contributions.
  virtual double externalEnergy(const Vectors& u, const TimeDomain& time) const;

public:
  //! \brief Prints a norm group to the log stream.
  //! \param[in] gNorm The norm values to print
  //! \param[in] rNorm Reference norms for the first norm group
  //! \param[in] prjName Projection name associated with this norm group
  virtual void printNormGroup (const Vector& gNorm, const Vector& rNorm,
                               const std::string& prjName) const;

private:
  std::vector<Material*> mVec;    //!< Material data
  RealArray              tVec;    //!< Beam thickness data
  PloadVec               myLoads; //!< Nodal point loads
};

#endif
