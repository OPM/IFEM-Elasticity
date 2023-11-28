// $Id$
//==============================================================================
//!
//! \file BeamProperty.h
//!
//! \date Apr 7 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class representing beam cross section properties.
//!
//==============================================================================

#ifndef _BEAM_PROPERTY_H
#define _BEAM_PROPERTY_H

class Vec3;
class RealFunc;
namespace tinyxml2 { class XMLElement; }


/*!
  \brief Class representing beam cross section properties.
*/

class BeamProperty
{
public:
  //! \brief Constructor initializing properties from an XML-element.
  explicit BeamProperty(const tinyxml2::XMLElement* prop = nullptr);
  //! \brief The destructor deallocates the property functions.
  ~BeamProperty();

  //! \brief Parses beam cross section properties from an XML-element.
  void parse(const tinyxml2::XMLElement* prop);
  //! \brief Parses circular cross section properties from an XML-element.
  static bool parsePipe(const tinyxml2::XMLElement* prop, double& A, double& I);
  //! \brief Parses massive box cross section properties from an XML-element.
  static bool parseBox(const tinyxml2::XMLElement* prop, double& A,
                       double& Iy, double& Iz);

  //! \brief Evaluates the beam properties at the specified point \a X.
  void eval(const Vec3& X, double L, double E, double G, double rho,
            bool hasGrav, bool hasMass,
            double& EA,   double& EI_y, double& EI_z,
            double& GI_t, double& Al_y, double& Al_z,
            double& rhoA, double& CG_y, double& CG_z,
            double& I_xx, double& I_yy, double& I_zz) const;
  //! \brief Evaluates the beam properties at the specified point \a X.
  void eval(const Vec3& X, double E, double G, double rho,
            double& EA, double& EI_y, double& EI_z, double& GI_t,
            double& rhoA, double& I_yy, double& I_zz) const;
  //! \brief Evaluates the line mass density at the specified point \a X.
  double evalRho(const Vec3& X, double rho) const;

  //! \brief Returns \e true unless all properties are constant.
  bool nonConstant() const
  {
    return (EAfunc || EIyfunc || EIzfunc || GItfunc ||
            rhofunc || Ixfunc || Iyfunc || Izfunc ||
            CGyfunc || CGzfunc);
  }

private:
  RealFunc* EAfunc;  //!< Axial stiffness
  RealFunc* EIyfunc; //!< Bending stiffness about local Y-axis
  RealFunc* EIzfunc; //!< Bending stiffness about local Z-axis
  RealFunc* GItfunc; //!< Torsional stiffness
  RealFunc* rhofunc; //!< Mass density (line mass)
  RealFunc* Ixfunc;  //!< Polar inertia
  RealFunc* Iyfunc;  //!< Inertia about local Y-axis
  RealFunc* Izfunc;  //!< Inertia about local Z-axis
  RealFunc* CGyfunc; //!< Mass center location along local Y-axis
  RealFunc* CGzfunc; //!< Mass center location along local Z-axis

public:
  double A;  //!< Cross section area
  double Ix; //!< Second area moment around local X-axis
  double Iy; //!< Second area moment around local Y-axis
  double Iz; //!< Second area moment around local Z-axis
  double It; //!< Torsional constant

  double Sy; //!< Shear centre offset in local Y-direction
  double Sz; //!< Shear centre offset in local Z-direction
  double Ky; //!< Shear reduction factor in local Y-direction
  double Kz; //!< Shear reduction factor in local Z-direction
};

#endif
