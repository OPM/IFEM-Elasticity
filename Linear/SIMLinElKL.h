// $Id$
//==============================================================================
//!
//! \file SIMLinElKL.h
//!
//! \date Sep 16 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based FE analysis of Kirchhoff-Love shells.
//!
//==============================================================================

#ifndef _SIM_LIN_EL_KL_H
#define _SIM_LIN_EL_KL_H

#include "SIMElasticity.h"
#include "SIM2D.h"


/*!
  \brief Driver class for isogeometric FEM analysis of Kirchhoff-Love shells.
*/

class SIMLinElKL : public SIMElasticity<SIM2D>
{
  /*!
    \brief Struct defining a nodal point load.
  */
  struct PointLoad
  {
    size_t patch; //!< Patch index [0,nPatch>
    int    inod;  //!< Local node number of the closest node
    char   ldof;  //!< Local DOF number of the load
    double xi[2]; //!< Parameters of the point (u,v)
    Vec3   X;     //!< Spatial coordinates of the point
    double pload; //!< Load magnitude
    // \brief Default constructor.
    PointLoad() : patch(0), inod(0), ldof(1) { xi[0] = xi[1] = pload = 0.0; }
  };

  typedef std::vector<PointLoad> PloadVec; //!< Point load container

public:
  //! \brief Default constructor.
  SIMLinElKL(bool shell = false);
  //! \brief Destructor.
  virtual ~SIMLinElKL();

  //! \brief Returns whether an analytical solution is available or not.
  virtual bool haveAnaSol() const;

protected:
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand();
  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Performs some pre-processing tasks on the FE model.
  virtual void preprocessA();
  //! \brief Performs some pre-processing tasks on the FE model.
  virtual bool preprocessB();

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd);
  //! \brief Initializes the body load properties for current patch.
  //! \param[in] patchInd 1-based patch index
  virtual bool initBodyLoad(size_t patchInd);
  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd);

  //! \brief Assembles the nodal point loads, if any.
  virtual bool assembleDiscreteTerms(const IntegrandBase*, const TimeDomain&);

  //! \brief Computes problem-dependent external energy contributions.
  virtual double externalEnergy(const Vectors& psol) const;

private:
  RealArray tVec;     //!< Shell thickness data
  PloadVec  myLoads;  //!< Nodal point loads
  int       aCode[3]; //!< Analytical BC codes (used by destructor)
};

#endif
