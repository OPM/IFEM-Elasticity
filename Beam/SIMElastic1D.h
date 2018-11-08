// $Id$
//==============================================================================
//!
//! \file SIMElastic1D.h
//!
//! \date Nov 10 2018
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for 1D elastic solution drivers.
//!
//==============================================================================

#ifndef _SIM_ELASTIC_1D_H
#define _SIM_ELASTIC_1D_H

#include "SIM1D.h"


/*!
  \brief Base class for 1D elastic solution drivers.
*/

class SIMElastic1D : public SIM1D
{
protected:
  //! \brief Default constructor.
  //! \param[in] n1 Dimension of the primary solution field
  explicit SIMElastic1D(unsigned char n1 = 1) : SIM1D(n1) {}
  //! \brief Empty destructor.
  virtual ~SIMElastic1D() {}

  //! \brief Finds the (closest) node or element of a specified point load.
  //! \param[in] ipt Load point identifier
  //! \param[in] patch 1-based index of the patch containing the load point
  //! \param xi Parameter value in range [0,1] of the load point.
  //! Converted to the parameter domain range of the found point on output,
  //! if \a onElement is \e true.
  //! \param[in] onElement If \e true, always find an element containing the
  //! specified point, even if it exactly matches a node (control point).
  //! \return 1-based patch-local index of the matching node, if positive.
  //! 1-based index of the element containing the point, if negative.
  //! If zero, an error condition occurred.
  int findLoadPoint(int ipt, int patch, double& xi, bool onElement);

  //! \brief Assembles consistent nodal forces due to an element point load.
  //! \param[in] patch 1-based patch index
  //! \param[in] u Parameter value of the loaded point
  //! \param[in] f Load magnitude
  //! \param[in] ldof Coordinate direction of the load
  bool assemblePoint(int patch, double u, double f, int ldof = 1);
};

#endif
