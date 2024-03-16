// $Id$
//==============================================================================
//!
//! \file ElasticBeam.h
//!
//! \date Aug 10 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class representing a nonlinear elastic beam.
//!
//==============================================================================

#ifndef _ELASTIC_BEAM_H
#define _ELASTIC_BEAM_H

#include "ElasticBase.h"

class BeamProperty;
class VecFunc;
namespace tinyxml2 { class XMLElement; }


/*!
  \brief Class representing the integrand of a nonlinear elastic beam.
*/

class ElasticBeam : public ElasticBase
{
public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of consequtive solution vectors to reside in core
  explicit ElasticBeam(unsigned short int n = 1);
  //! \brief The destructor deallocates the load and property functions.
  virtual ~ElasticBeam();

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Defines the (constant) beam stiffness.
  void setStiffness(double e, double g) { E = e; G = g; }
  //! \brief Defines the (constant) line mass.
  void setMass(double mass) { rho = mass; }

  //! \brief Assigns a load load function the the beam.
  void setBeamLoad(VecFunc* load, bool cpl = false);
  //! \brief Parses line load defintion from an XML-element.
  void parseBeamLoad(const tinyxml2::XMLElement* prop);
  //! \brief Parses beam cross section properties from an XML-element.
  static BeamProperty* parseProp(const tinyxml2::XMLElement* prop);
  //! \brief Defines the actual beam cross section properties to use.
  void setProperty(BeamProperty* prop) { myProp = prop; }

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const
  { return NO_DERIVATIVES | ELEMENT_CORNERS | NODAL_ROTATIONS; }

  using ElasticBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] iEl Global element number (1-based)
  virtual LocalIntegral* getLocalIntegral(size_t, size_t iEl, bool) const;

  using ElasticBase::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] fe Nodal and integration point data for current element
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC,
                           const FiniteElement& fe, const Vec3&, size_t,
                           LocalIntegral& elmInt);

  using ElasticBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  using ElasticBase::finalizeElement;
  //! \brief Finalizes the element matrices after the numerical integration.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] iGP Global integration point counter of first point in element
  virtual bool finalizeElement(LocalIntegral& elmInt,
                               const TimeDomain& time, size_t iGP);

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  virtual NormBase* getNormIntegrand(AnaSol* = nullptr) const;

  //! \brief Evaluates the displacement at current point.
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] eV Primary solution vector for current element
  virtual Vec3 displacement(const FiniteElement& fe, const Vector& eV) const;

protected:
  //! \brief Initializes the property function pointers.
  void initPropFunc();

  //! \brief Returns the local-to-global transformation matrix for an element.
  Matrix& getLocalAxes(LocalIntegral& elmInt) const;

  //! \brief Calculates the material stiffness matrix of the beam element.
  void getMaterialStiffness(Matrix& EK, double L,
                            double EA,  double GIt,
                            double EIy, double EIz,
                            double ALy, double ALz,
                            double Sy,  double Sz) const;
  //! \brief Calculates the geometric stiffness matrix of the beam element.
  void getGeometricStiffness(Matrix& EK, double N, double L,
                             double EIy, double EIz,
                             double ALy, double ALz, double ItoA,
                             double Sy,  double Sz) const;
  //! \brief Calculates the mass matrix of the beam element.
  void getMassMatrix(Matrix& EM, double rhoA, double Ixx,
                     double Iyy, double Izz, double L) const;

protected:
  // Material property parameters (constant)
  double E;   //!< Young's modulus
  double G;   //!< Shear modulus
  double rho; //!< Mass density

  BeamProperty* myProp;   //!< Geometric and constitutive beam properties
  VecFunc*      lineLoad; //!< Time and/or space-dependent line load
  VecFunc*      cplLoad;  //!< Time and/or space-dependent moment load

  //! If \e true, the element matrices are established in local axes
  //! and must be transformed to global axis directions before assembly
  bool inLocalAxes;
};


/*!
  \brief Class representing the integrand of elastic beam solution norms.
*/

class ElasticBeamNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The elastic beam problem to evaluate norms for
  //! \param[in] a The analytical displacement field
  ElasticBeamNorm(ElasticBeam& p, VecFunc* a) : NormBase(p), anasol(a) {}
  //! \brief Empty destructor.
  virtual ~ElasticBeamNorm() {}

  //! \brief Returns the number of reduced-order integration points.
  virtual int getReducedIntegration(int) const { return 0; }

  using NormBase::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC,
                           const FiniteElement&, const Vec3&, size_t,
                           LocalIntegral& elmInt);

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  virtual size_t getNoFields(int group = 0) const;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prfix Common prefix for all norm names
  virtual std::string getName(size_t, size_t j, const char* prfix) const;

private:
  VecFunc* anasol; //!< Analytical displacement field
};

#endif
