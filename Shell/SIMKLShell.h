// $Id$
//==============================================================================
//!
//! \file SIMKLShell.h
//!
//! \date Sep 16 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based FE analysis of Kirchhoff-Love shells.
//!
//==============================================================================

#ifndef _SIM_KL_SHELL_H
#define _SIM_KL_SHELL_H

#include "SIMElasticity.h"
#include "SIM2D.h"
#include "Interface.h"

class KirchhoffLove;

typedef std::pair<double,double> Doubles; //!< Convenience type


/*!
  \brief Driver class for isogeometric FEM analysis of Kirchhoff-Love shells.
*/

class SIMKLShell : public SIMElasticity<SIM2D>
{
protected:
  typedef std::pair<int,int> Ipair; //!< Convenience declaration

  /*!
    \brief Struct defining a nodal point load.
  */
  struct PointLoad
  {
    size_t      patch; //!< Patch index [0,nPatch>
    Ipair       ldof;  //!< Local node/DOF number of the closest node
    double      xi[2]; //!< Parameters of the point (u,v)
    Vec3        X;     //!< Spatial coordinates of the point
    ScalarFunc* p;     //!< Load magnitude
    //! \brief Default constructor.
    PointLoad() : patch(0), ldof(Ipair(0,1)), p(nullptr) { xi[0]=xi[1] = 0.0; }
  };

  typedef std::vector<PointLoad> PloadVec; //!< Point load container

  /*!
    \brief Struct defining a line load domain within a patch.
  */
  struct LLdomain
  {
    char    direction; //!< Line direction flag: 1=East-West, 2=North-South
    double  u;         //!< Location of the load line
    Doubles range;     //!< Parameter range of the load line
    //! \brief Default constructor.
    LLdomain() : direction(0) { u = range.first = range.second = 0.0; }
  };

public:
  //! \brief Default constructor.
  explicit SIMKLShell(const char* heading, bool isShell = true);
  //! \brief The destructor deletes the nodal point load functions.
  virtual ~SIMKLShell();

  //! \brief Returns that no analytical solution is available.
  virtual bool haveAnaSol() const { return false; }

  //! \brief Initializes the property containers of the model.
  virtual void clearProperties();

  //! \brief Prints a norm group to the log stream.
  //! \param[in] norm The norm values to print
  //! \param[in] rNorm Reference norms for the first norm group
  //! \param[in] prjName Projection name associated with this norm group
  virtual void printNormGroup(const Vector& norm, const Vector& rNorm,
                              const std::string& prjName) const;

protected:
  //! \brief Returns the actual integrand.
  virtual KirchhoffLove* getProblem(int version = 1);
  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand();
  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const tinyxml2::XMLElement* elem);

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

  //! \brief Returns the interface checker for line load terms in the integrand.
  virtual ASM::InterfaceChecker* getInterfaceChecker(size_t pidx) const;

  //! \brief Assembles the nodal point loads, if any.
  virtual bool assembleDiscreteTerms(const IntegrandBase*, const TimeDomain&);

  //! \brief Computes problem-dependent external energy contributions.
  virtual double externalEnergy(const Vectors& u, const TimeDomain& time) const;

  //! \brief Shifts global node numbers by a constant offset.
  virtual void shiftGlobalNums(int nshift, int);

  //! \brief Returns norm index of the integrated volume.
  virtual size_t getVolumeIndex() const { return 0; }

private:
  //! \brief Assembles consistent nodal forces due to an element point load.
  //! \param[in] patch 1-based patch index
  //! \param[in] u Parameter values of the loaded point
  //! \param[in] f Load magnitude
  //! \param[in] ldof Coordinate direction of the load
  bool assemblePoint(int patch, const double* u, double f, int ldof);

protected:
  RealArray tVec;     //!< Shell thickness data
  PloadVec  myLoads;  //!< Nodal/element point loads
  int       aCode[3]; //!< Analytical BC codes (used by destructor)
  LLdomain  lineLoad; //!< Domain definition of the line load
};


/*!
  \brief Checks for interface integration on the specified element boundaries.
*/

class ElmBorderChk
{
  double umin; //!< West border of load domain
  double umax; //!< East border of load domain
  double vmin; //!< South border of load domain
  double vmax; //!< North border of load domain

public:
  //! \brief Default constructor.
  ElmBorderChk();
  //! \brief Empty destructor.
  virtual ~ElmBorderChk() {}

  //! \brief Defines the parameter domain of a north-south line load.
  void setDomainNS(double u, const Doubles& rng);
  //! \brief Defines the parameter domain of an east-west line load.
  void setDomainEW(double v, const Doubles& rng);

  //! \brief Returns a status mask based on the element boundary parameters.
  //! \param[in] u0 Parameter value of the west element boundary
  //! \param[in] u1 Parameter value of the east element boundary
  //! \param[in] v0 Parameter value of the south element boundary
  //! \param[in] v1 Parameter value of the north element boundary
  short int maskBorder(double u0, double u1, double v0, double v1) const;
};


/*!
  \brief Checks for line load on the specified element boundaries.
*/

template<class T>
class LineLoadChecker : public T::InterfaceChecker, public ElmBorderChk
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  //! \param[in] pch Patch to check for line loads on
  explicit LineLoadChecker(const T& pch) : T::InterfaceChecker(pch) {}
  //! \brief Empty destructor.
  virtual ~LineLoadChecker() {}

  //! \brief Returns a status mask based on the element boundary parameters.
  //! \param[in] u0 Parameter value of the west element boundary
  //! \param[in] u1 Parameter value of the east element boundary
  //! \param[in] v0 Parameter value of the south element boundary
  //! \param[in] v1 Parameter value of the north element boundary
  virtual short int elmBorderMask(double u0, double u1,
                                  double v0, double v1, double, double) const
  {
    return this->maskBorder(u0,u1,v0,v1);
  }
};

#endif
