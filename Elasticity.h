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
class Tensor;
class SymmTensor;
class TractionFunc;
class FunctionBase;


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
  virtual Material* parseMatProp(const tinyxml2::XMLElement* elem,
                                 bool planeStrain = true);

  //! \brief Parses local coordinate system definition from a character string.
  bool parseLocalSystem(const char* cline);
  //! \brief Parses local coordinate system definition from an XML-element.
  bool parseLocalSystem(const tinyxml2::XMLElement* elem);

  //! \brief Parses a data section from an XML-element.
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(TractionFunc* tf) { tracFld = tf; }
  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(VecFunc* tf) { fluxFld = tf; }
  //! \brief Defines the body force field.
  void setBodyForce(VecFunc* bf) { bodyFld = bf; }
  //! \brief Defines the extraction function of the dual problem.
  void setDualRHS(FunctionBase* df) { dualRHS = df; }
  //! \brief Defines an extraction function for VCP.
  void addExtrFunction(FunctionBase* exf);
  //! \brief Returns the number of extraction functions.
  size_t numExtrFunction() const { return dualFld.size(); }

  //! \brief Defines the material properties.
  virtual void setMaterial(Material* mat);
  //! \brief Returns the current material object.
  Material* getMaterial() const { return material; }

  //! \brief Initializes time integration parameters.
  virtual bool init(const TimeDomain&) { return true; }

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
  //! \param[in] iEl Global element number (1-based)
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t iEl,
                                          bool neumann) const;

  //! \brief Defines the global integral for calculating reaction forces only.
  virtual void setSecondaryInt(GlobalIntegral* gq);
  //! \brief Returns the system quantity to be integrated by \a *this.
  virtual GlobalIntegral& getGlobalInt(GlobalIntegral* gq) const;

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const;

  using ElasticBase::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  using ElasticBase::evalSol;
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

  //! \brief Evaluates the finite element (FE) strain at an integration point.
  //! \param[out] s The FE strans values at current point
  //! \param[in] eV Element displacement vector
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  bool evalEps(Vector& s, const Vector& eV, const FiniteElement& fe,
               const Vec3& X) const;

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
  //! \param[in] X Cartesian coordinates of evaluation point
  //! \param[in] n Boundary normal vector at evaluation point
  //! \param[in] grd If \e true, evaluate the time-derivative of the traction
  Vec3 getTraction(const Vec3& X, const Vec3& n, bool grd = false) const;
  //! \brief Evaluates the body force field (if any) at specified point.
  //! \param[in] X Cartesian coordinates of evaluation point
  //! \param[in] grd If \e true, evaluate the time-derivative of the body force
  Vec3 getBodyforce(const Vec3& X, bool grd = false) const;
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

  //! \brief Returns the patch-wise extraction function field, if any.
  //! \param[in] ifield 1-based index of the field to return
  virtual Vector* getExtractionField(size_t ifield);

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution fields (optional)
  virtual NormBase* getNormIntegrand(AnaSol* asol = nullptr) const;

  //! \brief Returns a pointer to an Integrand for boundary force evaluation.
  //! \note The Integrand is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer returned.
  //! \param[in] x Reference point for torque calculation
  //! \param[in] asol Pointer to analytical solution fields (optional)
  virtual ForceBase* getForceIntegrand(const Vec3* x, AnaSol* asol) const;
  //! \brief Returns a pointer to an Integrand for nodal force evaluation.
  //! \note The Integrand is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer returned.
  virtual ForceBase* getForceIntegrand() const;

  //! \brief Returns norm index of the integrated volume.
  virtual size_t getVolumeIndex(bool withAnaSol) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld = 2) const;
  //! \brief Returns the name of a primary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField1Name(size_t i, const char* prefix) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField2Name(size_t i, const char* prefix) const;

  //! \brief Initializes the maximum stress values buffer.
  //! \param[in] nP Number of patches in the model, or 1 if global maximum
  //!
  //! \details This method allocates and/or initializes the internal array used
  //! to calculate the maximum stress values in the model. If the argument \a nP
  //! is equal to 1, only the global maxima over all patches in the model are
  //! determined. Otherwise, if it equals the number of patches in the model,
  //! the maximum values are computed for each patch separately.
  void initMaxVals(size_t nP = 1);

  //! \brief Returns a pointer to the max values for external update.
  std::vector<PointValues>* getMaxVals() const { return &maxVal; }

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
  //! \param sumLoad Total external load in each spatial direction
  //! \param[in] N Basis function values at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] grd If \e true, the gradient (time-derivative) is computed
  void formBodyForce(Vector& ES, RealArray& sumLoad, const Vector& N,
                     const Vec3& X, double detJW, bool grd = false) const;

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

  //! \brief Calculates the deformation gradient at current point.
  //! \param[in] eV Element solution vector
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] r Radial coordinate of current point
  //! \param[out] F Deformation gradient (or dUdX) at current point
  //! \param[in] gradOnly If \e true, only evaluate the gradient tensor dUdX
  bool formDefGradient(const Vector& eV, const Vector& N, const Matrix& dNdX,
                       double r, Tensor& F, bool gradOnly = false) const;

  //! \brief Evaluates the thermal strain at current integration point.
  virtual double getThermalStrain(const Vector&, const Vector&,
                                  const Vec3&) const { return 0.0; }

  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] s The FE stress values at current point
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol2(Vector& s, const Vectors& eV,
                        const FiniteElement& fe, const Vec3& X) const;

  //! \brief Performs pull-back of traction (interface for nonlinear problems).
  virtual bool pullBackTraction(Vec3&) const { return true; }

public:
  //! \brief Sets up the constitutive matrix at current point.
  //! \param[out] C \f$6\times6\f$-matrix (in 3D) or \f$3\times3\f$-matrix
  //! (in 2D), representing the constitutive tensor or its inverse
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] inverted If \e true, form the inverse constitutive tensor
  bool formCmat(Matrix& C, const FiniteElement& fe, const Vec3& X,
                bool inverted = false) const;

  //! \brief Returns \e true if this is an axial-symmetric problem.
  bool isAxiSymmetric() const { return axiSymmetry; }

  //! \brief Returns the tensile energy array (interface for fracture problems).
  virtual const RealArray* getTensileEnergy() const { return nullptr; }

  //! \brief Enable or disable max value calculation.
  void enableMaxValCalc(bool onOrOff) const { calcMaxVal = onOrOff; }

protected:
  // Physical properties
  Material*     material; //!< Material data and constitutive relation
  LocalSystem*  locSys;   //!< Local coordinate system for result output
  TractionFunc* tracFld;  //!< Pointer to implicit boundary traction field
  VecFunc*      fluxFld;  //!< Pointer to explicit boundary traction field
  VecFunc*      bodyFld;  //!< Pointer to body force field
  Vec3Vec*      pDirBuf;  //!< Principal stress directions buffer

  FunctionBase*              dualRHS; //!< Extraction function for dual RHS
  std::vector<FunctionBase*> dualFld; //!< Extraction functions for VCP

  mutable std::vector<PointValues> maxVal;  //!< Maximum result values
  mutable std::vector<Vec3Pair>    tracVal; //!< Traction field point values

  unsigned short int nDF; //!< Dimension on deformation gradient (2 or 3)
  bool       axiSymmetry; //!< If \e true, the problem is axi-symmetric
  double           gamma; //!< Numeric stabilization parameter

private:
  mutable bool calcMaxVal; //!< If \e true, max result values are calculated
  GlobalIntegral* myReacI; //!< Reaction-forces-only integral

public:
  static bool wantStrain;          //!< Option for output of strain, not stress
  static bool wantPrincipalStress; //!< Option for principal stress calculation
  static bool asolProject;         //!< Analytical solution projection flag
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
  //! \param[in] fld which field set to consider (2=all, 3=stress comps. only)
  explicit ElasticityNorm(Elasticity& p, STensorFunc* a = nullptr, int fld = 2);
  //! \brief Empty destructor.
  virtual ~ElasticityNorm() {}

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return true; }

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  using NormBase::evalBou;
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

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  virtual size_t getNoFields(int group) const;

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
  ElasticityForce(Elasticity& p, const Vec3* x, AnaSol* a = nullptr)
    : ForceBase(p), X0(x), anasol(a), nodal(false) {}
  //! \brief Constructor for global nodal force integration.
  //! \param[in] p The Elasticity problem to evaluate nodal forces for
  explicit ElasticityForce(Elasticity& p)
    : ForceBase(p), X0(nullptr), anasol(nullptr), nodal(true) {}

  //! \brief Empty destructor.
  virtual ~ElasticityForce() {}

  using ForceBase::getLocalIntegral;
  //! \brief Returns a local integral container for the element \a iEl.
  virtual LocalIntegral* getLocalIntegral(size_t, size_t iEl, bool) const;

  using ForceBase::evalBou;
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
