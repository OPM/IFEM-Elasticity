// $Id$
//==============================================================================
//!
//! \file SIMLinEl.h
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#ifndef _SIM_LIN_EL_H
#define _SIM_LIN_EL_H

#include "SIMElasticity.h"
#include "LinearElasticity.h"
#include "SIM2D.h"
#include "SIM3D.h"


/*!
  \brief Driver class for isogeometric FEM analysis of linear elastic problems.
*/

template<class Dim> class SIMLinEl : public SIMElasticity<Dim>
{
public:
  //! \brief Default constructor.
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  explicit SIMLinEl(bool checkRHS = false) : SIMElasticity<Dim>(checkRHS) {}
  //! \brief Empty destructor.
  virtual ~SIMLinEl() {}

protected:
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand()
  {
    if (!Dim::myProblem)
      Dim::myProblem = new LinearElasticity(Dim::dimension,
                                            Elastic::axiSymmetry,
                                            Elastic::GIpointsVTF);

    return dynamic_cast<Elasticity*>(Dim::myProblem);
  }

  //! \brief Parses a data section from an input file.
  //! \details This function allows for specialization of the template
  //! while still reusing as much code as possible.
  //! Only put dimension-specific code in here.
  virtual bool parseDimSpecific(char* keyWord, std::istream& is);

  //! \brief Parses a data section from an XML element.
  //! \details This function allows for specialization of the template
  //! while still reusing as much code as possible.
  //! Only put dimension-specific code in here.
  virtual bool parseDimSpecific(const TiXmlElement* elem);
};

typedef SIMLinEl<SIM2D> SIMLinEl2D; //!< 2D specific driver
typedef SIMLinEl<SIM3D> SIMLinEl3D; //!< 3D specific driver

//! \brief Template specialization - 2D specific input parsing.
template<> bool SIMLinEl2D::parseDimSpecific(char* keyWord, std::istream& is);
//! \brief Template specialization - 2D specific input parsing.
template<> bool SIMLinEl2D::parseDimSpecific(const TiXmlElement* elem);

//! \brief Template specialization - 3D specific input parsing.
template<> bool SIMLinEl3D::parseDimSpecific(char* keyWord, std::istream& is);
//! \brief Template specialization - 3D specific input parsing.
template<> bool SIMLinEl3D::parseDimSpecific(const TiXmlElement* elem);

#endif
