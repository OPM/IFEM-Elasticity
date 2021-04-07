// $Id$
//==============================================================================
//!
//! \file SIMElasticBar.h
//!
//! \date Aug 11 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based FE analysis of elastic bars & beams.
//!
//==============================================================================

#ifndef _SIM_ELASTIC_BAR_H
#define _SIM_ELASTIC_BAR_H

#include "SIMElastic1D.h"

class ElasticBar;
class ElasticBeam;
class BeamProperty;
class Tensor;


/*!
  \brief Driver class for isogeometric FEM analysis of elastic bars and beams.
*/

class SIMElasticBar : public SIMElastic1D
{
public:
  //! \brief Default constructor.
  //! \param[in] n Number of consequtive solution vectors in core
  explicit SIMElasticBar(unsigned char n = 1);
  //! \brief The destructor deletes the nodal point load functions.
  virtual ~SIMElasticBar();

  //! \brief Prints out problem-specific data to the log stream.
  virtual void printProblem() const;

  //! \brief Creates the computational FEM model from the spline patches.
  //! \details Reimplemented to account for twist angle in beam problems.
  virtual bool createFEMmodel(char);

  //! \brief Updates the nodal rotations for problems with rotational DOFs.
  //! \param[in] incSol Incremental solution to update the rotations with
  //! \param[in] alpha Scaling factor for the incremental solution.
  //! If 0.0, reinitialize the rotations from unity
  virtual bool updateRotations(const Vector& incSol, double alpha);

  //! \brief Returns the current rotation tensor for the specified global node.
  Tensor getNodeRotation(int inod) const;

  //! \brief Computes the total external load of current load/time step.
  //! \param[out] extloa Sum external load in each direction
  //! \param[in] time Parameters for nonlinear simulations
  virtual bool getExtLoad(RealArray& extloa, const TimeDomain& time) const;

protected:
  //! \brief Returns the actual beam problem integrand.
  virtual ElasticBar* getBarIntegrand(const std::string&);
  //! \brief Returns the actual beam problem integrand.
  virtual ElasticBeam* getBeamIntegrand(const std::string&);

  //! \brief Returns a list of prioritized XML-tags.
  virtual const char** getPrioritizedTags() const;

  using SIM1D::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Parses the twist angle description along the curve.
  bool parseTwist(const TiXmlElement* elem);

  //! \brief Preprocessing performed before the FEM model generation.
  virtual void preprocessA();
  //! \brief Preprocessing performed after the FEM model generation.
  virtual bool preprocessB();

  //! \brief Initializes beam properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd);
  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd);

  //! \brief Renumbers all global nodes number if the model.
  //! \param[in] nodeMap Mapping from old to new node number
  virtual bool renumberNodes(const std::map<int,int>& nodeMap);

  //! \brief Assembles the nodal point loads, if any.
  virtual bool assembleDiscreteTerms(const IntegrandBase* itg,
                                     const TimeDomain& time);

private:
  //! \brief Struct defining a nodal point load.
  struct PointLoad
  {
    int         inod; //!< Node or patch index
    int         ldof; //!< Local DOF number
    double      xi;   //!< Parameter of the point
    ScalarFunc* p;    //!< Load magnitude
    //! \brief Default constructor.
    PointLoad(int n = 0) : inod(n), ldof(0), xi(-1.0), p(nullptr) {}
  };

  std::vector<PointLoad>     myLoads; //!< Nodal/element point loads
  std::vector<BeamProperty*> myBCSec; //!< Beam cross section properties

  mutable bool printed; //!< If \e true, the problem definition as been printed
  char         lcStiff; //!< Flag for inclusion of load correction stiffness

protected:
  unsigned char nsv; //!< Number of consequtive solution vectors in core

  RealFunc* twist; //!< Twist angle along the beam
  Vec3      XZp;   //!< Local Z-direction vector
};

#endif
