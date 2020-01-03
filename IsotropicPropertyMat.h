// $Id$
//==============================================================================
//!
//! \file IsotropicPropertyMat.h
//!
//! \date Sep 13 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Isotropic linear elastic material model defined as a table.
//!
//==============================================================================

#ifndef _ISOTROPIC_PROPERTY_MAT_H
#define _ISOTROPIC_PROPERTY_MAT_H

#include "Function.h"
#include "LinIsotropic.h"
#include "TextureProperties.h"
#include <array>
#include <map>


/*!
  \brief Class representing an isotropic material model with a texture.
*/

class IsotropicPropertyMat : public LinIsotropic, public TextureProperties
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  //! \param[in] ps If \e true, assume plane stress in 2D
  //! \param[in] ax If \e true, assume 3D axi-symmetric material
  IsotropicPropertyMat(bool ps, bool ax) : LinIsotropic(ps,ax) {}
  //! \brief Empty destructor.
  virtual ~IsotropicPropertyMat() = default;

  //! \brief Parses material parameters from an XML element.
  void parse(const TiXmlElement* elem) override;

  //! \brief Prints out material parameters to the log stream.
  void printLog() const override;

private:
  class PropertyFunc : public RealFunc {
  public:
    PropertyFunc(const std::string& prop,
                 const TextureProperties& props)
      : m_prop(prop), m_props(props)
    {}

    double evaluate(const Vec3& X) const override
    {
        double val;
        m_props.getProperty(m_prop, X, val);
        return val;
    }

  protected:
    std::string m_prop;
    const TextureProperties& m_props;
  };
};

#endif
