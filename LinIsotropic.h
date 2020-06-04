// $Id$
//==============================================================================
//!
//! \file LinIsotropic.h
//!
//! \date Mar 01 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Isotropic linear elastic material model.
//!
//==============================================================================

#ifndef _LIN_ISOTROPIC_H
#define _LIN_ISOTROPIC_H

#include "MaterialBase.h"
#include "Field.h"

class ScalarFunc;
class RealFunc;


/*!
  \brief Class representing an isotropic linear elastic material model.
*/

class LinIsotropic : public Material
{
public:
  //! \brief Default constructor.
  //! \param[in] ps If \e true, assume plane stress in 2D
  //! \param[in] ax If \e true, assume 3D axi-symmetric material
  LinIsotropic(bool ps = false, bool ax = false);
  //! \brief Constructor initializing the material parameters.
  //! \param[in] E Young's modulus (constant)
  //! \param[in] v Poisson's ratio
  //! \param[in] densty Mass density
  //! \param[in] ps If \e true, assume plane stress in 2D
  //! \param[in] ax If \e true, assume 3D axi-symmetric material
  LinIsotropic(double E, double v = 0.0, double densty = 0.0,
               bool ps = false, bool ax = false)
    : Efunc(nullptr), Efield(nullptr), Emod(E), nuFunc(nullptr), nu(v),
      rhoFunc(nullptr), rho(densty),
      Cpfunc(nullptr), heatcapacity(0.0), Afunc(nullptr), alpha(0.0),
      condFunc(nullptr), conductivity(0.0), planeStress(ps), axiSymmetry(ax) {}
  //! \brief Constructor initializing the material parameters.
  //! \param[in] E Young's modulus (spatial function)
  //! \param[in] v Poisson's ratio
  //! \param[in] density Mass density
  //! \param[in] ps If \e true, assume plane stress in 2D
  //! \param[in] ax If \e true, assume 3D axi-symmetric material
  LinIsotropic(RealFunc* E, double v = 0.0, double density = 0.0,
               bool ps = false, bool ax = false);
  //! \brief Constructor initializing the material parameters.
  //! \param[in] E Young's modulus (spatial field)
  //! \param[in] v Poisson's ratio
  //! \param[in] density Mass density
  //! \param[in] ps If \e true, assume plane stress in 2D
  //! \param[in] ax If \e true, assume 3D axi-symmetric material
  LinIsotropic(Field* E, double v = 0.0, double density = 0.0,
               bool ps = false, bool ax = false);
  //! \brief The destructor deletes the stiffness function, if defined.
  virtual ~LinIsotropic();

  //! \brief Parses material parementers from an XML element.
  virtual void parse(const TiXmlElement* elem);

  //! \brief Prints out material parameters to the log stream.
  virtual void printLog() const;

  //! \brief Returns \e false if plane stress in 2D.
  virtual bool isPlaneStrain() const { return !planeStress; }

  //! \brief Evaluates the stiffness at current point.
  virtual double getStiffness(const Vec3& X) const;
  //! \brief Evaluates the plate stiffness parameter at current point.
  virtual double getPlateStiffness(const Vec3& X, double t) const;
  //! \brief Evaluates the mass density at current point.
  virtual double getMassDensity(const Vec3&) const;
  //! \brief Evaluates the heat capacity for given temperature.
  virtual double getHeatCapacity(double T) const;
  //! \brief Evaluates the thermal conductivity for given temperature.
  virtual double getThermalConductivity(double T) const;
  //! \brief Evaluates the thermal expansion coefficient for given temperature.
  virtual double getThermalExpansion(double T) const;

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
  virtual bool evaluate(Matrix& C, SymmTensor& sigma, double& U,
                        const FiniteElement& fe, const Vec3& X,
                        const Tensor&, const SymmTensor& eps,
                        char iop = 1, const TimeDomain* = nullptr,
                        const Tensor* = nullptr) const;

  //! \brief Evaluates the Lame-parameters at an integration point.
  //! \param[out] lambda Lame's first parameter
  //! \param[out] mu Lame's second parameter (shear modulus)
  //! \param[in] fe Finite element quantities at current point
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evaluate(double& lambda, double& mu,
                        const FiniteElement& fe, const Vec3& X) const;

  //! \brief Returns the function, if any, describing the stiffness variation.
  const RealFunc* getEfunc() const { return Efunc; }
  //! \brief Returns the field, if any, describing the stiffness variation.
  const Field* getEfield() const { return Efield; }

protected:
  // Material properties
  RealFunc* Efunc;      //!< Young's modulus (spatial function)
  Field*    Efield;     //!< Young's modulus (spatial field)
  double    Emod;       //!< Young's modulus (constant)
  RealFunc* nuFunc;     //!< Poisson's ratio (spatial function)
  double    nu;         //!< Poisson's ratio (constant)
  RealFunc* rhoFunc;    //!< Mass density (spatial function)
  double    rho;        //!< Mass density (constant)
  ScalarFunc* Cpfunc;   //!< Specific heat capacity function
  double heatcapacity;  //!< Specific heat capacity (constant)
  ScalarFunc* Afunc;    //!< Thermal expansion coefficient function
  double alpha;         //!< Thermal expansion coefficient (constant)
  ScalarFunc* condFunc; //!< Thermal conductivity function
  double conductivity;  //!< Thermal conductivity (constant)
  bool   planeStress;   //!< Plane stress/strain option for 2D problems
  bool   axiSymmetry;   //!< Axi-symmetric option
};

#endif
