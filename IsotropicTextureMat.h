// $Id$
//==============================================================================
//!
//! \file IsotropicTextureMat.h
//!
//! \date Sep 13 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Isotropic linear elastic material model defined through a texture.
//!
//==============================================================================

#ifndef _ISOTROPIC_TEXTURE_MAT_H
#define _ISOTROPIC_TEXTURE_MAT_H

#include "LinIsotropic.h"
#include <array>
#include <map>


/*!
  \brief Class representing an isotropic material model with a texture.
*/

class IsotropicTextureMat : public LinIsotropic
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  //! \param[in] ps If \e true, assume plane stress in 2D
  //! \param[in] ax If \e true, assume 3D axi-symmetric material
  IsotropicTextureMat(bool ps, bool ax) : LinIsotropic(ps,ax) {}
  //! \brief Empty destructor.
  virtual ~IsotropicTextureMat() = default;

  //! \brief Parses material parameters from an XML element.
  void parse(const TiXmlElement* elem) override;

  //! \brief Prints out material parameters to the log stream.
  void printLog() const override;

  //! \brief Evaluates the constitutive relation at an integration point.
  //! \param[out] C Constitutive matrix at current point
  //! \param[out] sigma Stress tensor at current point
  //! \param[out] U Strain energy density at current point
  //! \param[in] fe Finite element quantities at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] eps Strain tensor at current point
  //! \param[in] iop Calculation option;
  //!  -1 : Calculate the inverse constitutive matrix only,
  //!   0 : Calculate the constitutive matrix only,
  //!   1 : Calculate Cauchy stresses and the constitutive matrix,
  //!   3 : Calculate the strain energy density only.
  bool evaluate(Matrix& C, SymmTensor& sigma, double& U,
                const FiniteElement& fe, const Vec3& X,
                const Tensor&, const SymmTensor& eps,
                char iop = 1, const TimeDomain* = nullptr,
                const Tensor* = nullptr) const override;

  //! \brief Evaluates the Lame-parameters at an integration point.
  //! \param[out] lambda Lame's first parameter
  //! \param[out] mu Lame's second parameter (shear modulus)
  //! \param[in] fe Finite element quantities at current point
  //! \param[in] X Cartesian coordinates of current point
  bool evaluate(double& lambda, double& mu,
                const FiniteElement& fe, const Vec3& X) const override;

private:
  //! \brief Locates the appropriate material as indicated by texture.
  const LinIsotropic* findMaterial(const FiniteElement& fe) const;

protected:
  typedef std::pair<double,double> Doubles; //!< Convenience type
  typedef std::array<double,4>     rgba;    //!< RGB color code

  //! Material for different texture regions
  std::map< Doubles,LinIsotropic > materials;
  //! Raw image texture information describing the spatial material variation
  std::vector< std::vector<rgba> > textureData;
};

#endif
