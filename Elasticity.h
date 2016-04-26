// $Id$
//==============================================================================
//!
//! \file Elasticity.h
//!
//! \date Nov 12 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for linear and nonlinear elasticity problems.
//!
//==============================================================================

#ifndef _ELASTICITY_H
#define _ELASTICITY_H

#include "ElasticBase.h"

class LocalSystem;
class Material;
class ElmNorm;
class ElmMats;
class TiXmlElement;


/*!
  \brief Base class representing the integrand of elasticity problems.
  \details Implements common features for linear and nonlinear elasticity
  problems. This class is used for continuum problems only (2D and 3D domains
  with the same number of unknowns per node as the number of space dimensions.
  Note that the \a evalInt method is not implemented by this class.
  Thus, it is regarded as an abstract base class with a protected constructor.
*/

class Elasticity : public ElasticBase
{
protected:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] ax \e If \e true, and axisymmetric 3D formulation is assumed
  Elasticity(unsigned short int n = 3, bool ax = false);

public:
  //! \brief The destructor frees the dynamically allocated data objects.
  virtual ~Elasticity();

  //! \brief Parses material properties from a character string.
  virtual Material* parseMatProp(char* cline, bool planeStrain = true);
  //! \brief Parses material properties from an XML-element.
  virtual Material* parseMatProp(const TiXmlElement* elem,
                                 bool planeStrain = true);

  //! \brief Parses local coordinate system definition from a character string.
  bool parseLocalSystem(const char* cline);
  //! \brief Parses local coordinate system definition from an XML-element.
  bool parseLocalSystem(const TiXmlElement* elem);

  //! \brief Parses a data section from an XML-element.
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(TractionFunc* tf) { tracFld = tf; }
  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(VecFunc* tf) { fluxFld = tf; }
  //! \brief Defines the body force field.
  void setBodyForce(VecFunc* bf) { bodyFld = bf; }

  //! \brief Defines the material properties.
  virtual void setMaterial(Material* mat) { material = mat; }
  //! \brief Defines the local coordinate system for stress output.
  void setLocalSystem(LocalSystem* cs) { locSys = cs; }

  using ElasticBase::initIntegration;
  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  //! \param[in] nBp Total number of boundary integration points
  virtual void initIntegration(size_t nGp, size_t nBp);

  //! \brief Initializes the integrand for a new result point loop.
  //! \param[in] lambda Load parameter
  //! \param[in] prinDirs If \e true, compute/store principal directions
  virtual void initResultPoints(double lambda, bool prinDirs);

  using ElasticBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return true; }

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  using ElasticBase::evalSol;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const FiniteElement& fe, const Vec3& X,
                       const std::vector<int>& MNPC) const;

  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] s The FE stress values at current point
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] toLocal If \e true, transform to local coordinates (if defined)
  //! \param[out] pdir Directions of the principal stresses (optional)
  virtual bool evalSol(Vector& s, const Vectors& eV, const FiniteElement& fe,
                       const Vec3& X, bool toLocal = false,
                       Vec3* pdir = nullptr) const;

  //! \brief Evaluates the analytical solution at an integration point.
  //! \param[out] s The analytical stress values at current point
  //! \param[in] asol The analytical solution field
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s, const STensorFunc& asol, const Vec3& X) const;

  //! \brief Evaluates the primary solution at a result point.
  //! \param[in] eV Element solution vector
  //! \param[in] N Basis function values at current point
  //! \return Primary solution vector at current point
  Vec3 evalSol(const Vector& eV, const Vector& N) const;

  //! \brief Returns an evaluated principal direction vector field for plotting.
  //! \param[out] pdir Principal direction vectors at each result point
  //! \param[in] nPt Number of result points
  //! \param[in] idx 1-based index of which direction vector to return
  virtual bool getPrincipalDir(Matrix& pdir, size_t nPt, size_t idx) const;

  //! \brief Evaluates the boundary traction field (if any) at specified point.
  Vec3 getTraction(const Vec3& X, const Vec3& n) const;
  //! \brief Evaluates the body force field (if any) at specified point.
  virtual Vec3 getBodyforce(const Vec3& X) const;
  //! \brief Returns whether an external load is defined.
  virtual bool haveLoads() const;

  //! \brief Writes the surface tractions for a given time step to VTF-file.
  //! \param vtf The VTF-file object to receive the tractions
  //! \param[in] iStep Load/time step identifier
  //! \param geoBlk Running geometry block counter
  //! \param nBlock Running result block counter
  virtual bool writeGlvT(VTF* vtf, int iStep, int& geoBlk, int& nBlock) const;

  //! \brief Returns whether there are any traction values to write to VTF.
  virtual bool hasTractionValues() const { return !tracVal.empty(); }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution fields (optional)
  virtual NormBase* getNormIntegrand(AnaSol* asol = 0) const;

  //! \brief Returns a pointer to an Integrand for boundary force evaluation.
  //! \note The Integrand is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer returned.
  //! \param[in] x Reference point for torque calculation
  //! \param[in] asol Pointer to analytical solution fields (optional)
  virtual ForceBase* getForceIntegrand(const Vec3* x, AnaSol* asol = 0) const;
  //! \brief Returns a pointer to an Integrand for nodal force evaluation.
  //! \note The Integrand is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer returned.
  virtual ForceBase* getForceIntegrand() const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld = 2) const;
  //! \brief Returns the name of a primary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField1Name(size_t i, const char* prefix = 0) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField2Name(size_t i, const char* prefix = 0) const;

  typedef std::pair<Vec3,double> PointValue; //!< Convenience type

  //! \brief Returns a pointer to the max values for external update.
  std::vector<PointValue>* getMaxVals() const { return &maxVal; }

  //! \brief Prints out the maximum secondary solution values to the log stream.
  //! \param[in] precision Number of digits after the decimal point
  //! \param[in] comp Which component to print (0 means all)
  void printMaxVals(std::streamsize precision, size_t comp = 0) const;

protected:
  //! \brief Calculates some kinematic quantities at current point.
  //! \param[in] eV Element solution vector
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] r Radial coordinate of current point
  //! \param[out] B The strain-displacement matrix
  //! \param[out] eps Strain tensor at current point
  virtual bool kinematics(const Vector& eV,
			  const Vector& N, const Matrix& dNdX, double r,
			  Matrix& B, Tensor&, SymmTensor& eps) const;

  //! \brief Calculates integration point geometric stiffness contributions.
  //! \param EM Element matrix to receive the stiffness contributions
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] r Radial coordinate of current point
  //! \param[in] sigma Stress tensor at current point
  //! \param[in] detJW Jacobian determinant times integration point weight
  void formKG(Matrix& EM, const Vector& N, const Matrix& dNdX,
	      double r, const Tensor& sigma, double detJW) const;

  //! \brief Calculates integration point mass matrix contributions.
  //! \param EM Element matrix to receive the mass contributions
  //! \param[in] N Basis function values at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] detJW Jacobian determinant times integration point weight
  void formMassMatrix(Matrix& EM, const Vector& N,
		      const Vec3& X, double detJW) const;

  //! \brief Calculates integration point body force vector contributions.
  //! \param ES Element vector to receive the body force contributions
  //! \param[in] N Basis function values at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] detJW Jacobian determinant times integration point weight
  void formBodyForce(Vector& ES, const Vector& N,
		     const Vec3& X, double detJW) const;

  //! \brief Calculates the strain-displacement matrix.
  //! \param[in] Bmat The strain-displacement matrix
  //! \param[in] dNdX Basis function gradients at current point
  bool formBmatrix(Matrix& Bmat, const Matrix& dNdX) const;
  //! \brief Calculates the axi-symmetric strain-displacement matrix.
  //! \param[in] Bmat The strain-displacement matrix
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] r Radial coordinate of current point
  bool formBmatrix(Matrix& Bmat, const Vector& N, const Matrix& dNdX,
		   double r) const;

  //! \brief Evaluates the thermal strain at current integration point.
  virtual double getThermalStrain(const Vector&, const Vector&,
                                  const Vec3&) const { return 0.0; }

  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] s The FE stress values at current point
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  bool evalSol2(Vector& s, const Vectors& eV,
                const FiniteElement& fe, const Vec3& X) const;

  //! \brief Performs pull-back of traction (interface for nonlinear problems).
  virtual bool pullBackTraction(Vec3&) const { return true; }

public:
  //! \brief Sets up the inverse constitutive matrix at current point.
  //! \param[out] Cinv \f$6\times6\f$-matrix (in 3D) or \f$3\times3\f$-matrix
  //! (in 2D), representing the inverse constitutive tensor
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  bool formCinverse(Matrix& Cinv, const FiniteElement& fe, const Vec3& X) const;

  //! \brief Returns \e true if this is an axial-symmetric problem.
  bool isAxiSymmetric() const { return axiSymmetry; }

  //! \brief Returns the tensile energy array (interface for fracture problems).
  virtual const RealArray* getTensileEnergy() const { return nullptr; }

protected:
  // Physical properties
  Material*     material; //!< Material data and constitutive relation
  LocalSystem*  locSys;   //!< Local coordinate system for result output
  TractionFunc* tracFld;  //!< Pointer to implicit boundary traction field
  VecFunc*      fluxFld;  //!< Pointer to explicit boundary traction field
  VecFunc*      bodyFld;  //!< Pointer to body force field
  Vec3Vec*      pDirBuf;  //!< Principal stress directions buffer

  mutable std::vector<PointValue> maxVal;  //!< Maximum result values
  mutable std::vector<Vec3Pair>   tracVal; //!< Traction field point values

  unsigned short int nDF; //!< Dimension on deformation gradient (2 or 3)
  bool       axiSymmetry; //!< \e true if the problem is axi-symmetric
  double           gamma; //!< Numeric stabilization parameter

public:
  static bool wantPrincipalStress; //!< Option for principal stress calculation
};


/*!
  \brief Class representing the integrand of elasticity energy norms.
*/

class ElasticityNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The linear elasticity problem to evaluate norms for
  //! \param[in] a The analytical stress field (optional)
  ElasticityNorm(Elasticity& p, STensorFunc* a = 0);
  //! \brief Empty destructor.
  virtual ~ElasticityNorm() {}

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return true; }

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  using NormBase::finalizeElement;
  //! \brief Finalizes the element norms after the numerical integration.
  //! \param elmInt The local integral object to receive the contributions
  //!
  //! \details This method is used to compute effectivity indices.
  virtual bool finalizeElement(LocalIntegral& elmInt);

  //! \brief Adds external energy terms to relevant norms.
  //! \param gNorm Global norm quantities
  //! \param[in] energy Global external energy
  virtual void addBoundaryTerms(Vectors& gNorm, double energy) const;

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  virtual size_t getNoFields(int group = 0) const;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  virtual std::string getName(size_t i, size_t j, const char* prefix) const;

  //! \brief Returns whether a norm quantity stores element contributions.
  virtual bool hasElementContributions(size_t i, size_t j) const;

private:
  STensorFunc* anasol; //!< Analytical stress field
};


/*!
  \brief Class representing the integrand for computing boundary forces
*/

class ElasticityForce : public ForceBase
{
public:
  //! \brief Constructor for global force resultant integration.
  //! \param[in] p The Elasticity problem to evaluate forces for
  //! \param[in] x Reference point for torque calculation
  //! \param[in] a The analytical velocity and pressure fields
  ElasticityForce(Elasticity& p, const Vec3* x, AnaSol* a = 0)
    : ForceBase(p), X0(x), anasol(a), nodal(false) {}
  //! \brief Constructor for global nodal force integration.
  //! \param[in] p The Elasticity problem to evaluate nodal forces for
  ElasticityForce(Elasticity& p) : ForceBase(p), X0(0), anasol(0), nodal(true) {}

  //! \brief Empty destructor.
  virtual ~ElasticityForce() {}

  //! \brief Returns a local integral container for the element \a iEl.
  virtual LocalIntegral* getLocalIntegral(size_t, size_t iEl, bool) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  //! \brief Returns the number of force components.
  virtual size_t getNoComps() const;

private:
  //! \brief Evaluates the integrand for global force resultants.
  bool evalForce(ElmNorm& pnorm, const Vec3& th,
                 const Vec3& X, const Vec3& normal,
                 double detJW, size_t nsd) const;

  //! \brief Evaluates the integrand for nodal forces.
  bool evalForce(ElmMats& elmat, const Vec3& th, const FiniteElement& fe) const;

  //! \brief Evaluates the integrand for Robin type nodal conditions.
  bool evalRobin(ElmMats& elmat, const Vec3& th, const Vec3& U,
                 const FiniteElement& fe, double alpha = 1.0) const;

protected:
  const Vec3* X0; //!< Reference point for torque computations
  AnaSol* anasol; //!< Analytical solution fields
  bool     nodal; //!< If \e true, we are integrating nodal forces
};

#endif
