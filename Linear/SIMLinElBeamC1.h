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

#include "SIMElastic1D.h"

class Material;


/*!
  \brief Driver class for isogeometric FEM analysis of C1-continuous beams.
*/

class SIMLinElBeamC1 : public SIMElastic1D
{
public:
  //! \brief Default constructor.
  explicit SIMLinElBeamC1(const char* h = nullptr) { if (h) myHeading = h; }
  //! \brief The destructor deletes the material data objects.
  virtual ~SIMLinElBeamC1();

  //! \brief Initializes the property containers of the model.
  virtual void clearProperties();

protected:
  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Parses or generates app-specific explicit knots for refinement.
  //! \param[in] elem The XML element to parse
  //! \param[out] xi Explicit knots in range [0,1]
  virtual bool parseXi(const tinyxml2::XMLElement* elem, RealArray& xi) const;

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

  //! \brief Shifts global node numbers by a constant offset.
  virtual void shiftGlobalNums(int nshift, int);

public:
  //! \brief Prints a norm group to the log stream.
  //! \param[in] gNorm The norm values to print
  //! \param[in] rNorm Reference norms for the first norm group
  //! \param[in] prjName Projection name associated with this norm group
  virtual void printNormGroup(const Vector& gNorm, const Vector& rNorm,
                              const std::string& prjName) const;

  //! \brief Returns whether an analytical solution is available or not.
  virtual bool haveAnaSol() const;

private:
  //! \brief Struct defining a point load.
  struct PointLoad
  {
    size_t patch; //!< Patch index [0,nPatch>
    int    inod;  //!< Local node/element number
    double xi;    //!< Parameter of the point
    double pload; //!< Load magnitude
    // \brief Default constructor.
    PointLoad() : patch(1), inod(0), xi(-1.0), pload(0.0) {}
  };

  std::vector<PointLoad> myLoads; //!< Nodal/element point loads
  std::vector<Material*> mVec;    //!< Material data
  RealArray              tVec;    //!< Beam thickness data
};

#endif
