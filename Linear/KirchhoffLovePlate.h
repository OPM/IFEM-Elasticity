// $Id$
//==============================================================================
//!
//! \file KirchhoffLovePlate.h
//!
//! \date Sep 13 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class for linear Kirchhoff-Love thin plate problems.
//!
//==============================================================================

#ifndef _KIRCHHOFF_LOVE_PLATE_H
#define _KIRCHHOFF_LOVE_PLATE_H

#include "KirchhoffLove.h"

class VecFunc;
class STensorFunc;


/*!
  \brief Class representing the integrand of thin plate problems.

  \details The formulation is based on Kirchhoff-Love plate theory
  and therefore requires second-derivatives of the basis functions.

  See the document doc/Integrands/KirchhoffLove.pdf
  for the theoretical foundation of the integrand implemented in
  the evalInt() and evalBou() methods.
*/

class KirchhoffLovePlate : public KirchhoffLove
{
public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions (1=beam, 2=plate)
  //! \param[in] v Integrand version (1: B-matrix, 2: Tensor form)
  explicit KirchhoffLovePlate(unsigned short int n = 2, short int v = 1);
  //! \brief Empty destructor,
  virtual ~KirchhoffLovePlate() {}

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Defines the Neumann order that is the subject of integration.
  virtual void setNeumannOrder(char flag) { nOrder = flag; }

  //! \brief Returns the integrand version flag.
  short int getVersion() const { return version; }

  using KirchhoffLove::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  using KirchhoffLove::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  using KirchhoffLove::evalSol;
  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] s The FE stress resultant values at current point
  //! \param[in] eV Element solution vector
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] toLocal If \e true, transform to local coordinates (if defined)
  virtual bool evalSol(Vector& s, const Vectors& eV,
                       const FiniteElement& fe, const Vec3& X,
                       bool toLocal = false) const;

  //! \brief Returns the plate stiffness parameter at the specified point \a X.
  double getStiffness(const Vec3& X) const;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution fields (optional)
  virtual NormBase* getNormIntegrand(AnaSol* asol = nullptr) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld = 2) const;
  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  virtual std::string getField1Name(size_t, const char* prefix) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField2Name(size_t i, const char* prefix) const;

protected:
  //! \brief Calculates the strain-displacement matrix \b B at current point.
  //! \param[out] Bmat The strain-displacement matrix
  //! \param[in] d2NdX2 Basis function 2nd derivatives at current point
  bool formBmatrix(Matrix& Bmat, const Matrix3D& d2NdX2) const;

  //! \brief Evaluates the stiffness matrix integrand, version 1.
  //! \param EK The element stiffness matrix to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalK1(Matrix& EK, const FiniteElement& fe, const Vec3& X) const;
  //! \brief Evaluates the stiffness matrix integrand, version 2.
  //! \param EK The element stiffness matrix to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalK2(Matrix& EK, const FiniteElement& fe, const Vec3& X) const;

public:
  //! \brief Sets up the constitutive matrix at current point.
  //! \param[out] C \f$3\times3\f$-matrix, representing the constitutive tensor
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] invers If \e true, the inverse matrix is establised instead
  bool formCmatrix(Matrix& C, const FiniteElement& fe, const Vec3& X,
                   bool invers = false) const;

protected:
  short int version; //!< Integrand version flag
  short int nOrder;  //!< Neumann order flag, 1 = moments, 2 = shear forces
};


/*!
  \brief Class representing the integrand of Kirchhoff-Love energy norms.
*/

class KirchhoffLovePlateNorm : public NormBase
{
public:
  //! \brief The constructor initializes its data members.
  //! \param[in] p The thin plate problem to evaluate norms for
  //! \param[in] a The analytical stress resultant field (optional)
  KirchhoffLovePlateNorm(KirchhoffLovePlate& p, STensorFunc* a = nullptr);
  //! \brief This constructor also initializes its data members.
  //! \param[in] p The thin plate problem to evaluate norms for
  //! \param[in] a The analytical 2nd derivatives of the displacement field
  KirchhoffLovePlateNorm(KirchhoffLovePlate& p, VecFunc* a);
  //! \brief Empty destructor.
  virtual ~KirchhoffLovePlateNorm() {}

  //! \brief Defines the Neumann order that is the subject of integration.
  virtual void setNeumannOrder(char flag) { nOrder = flag; }

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

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const;

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return true; }

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  virtual size_t getNoFields(int group = 0) const;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  virtual std::string getName(size_t i, size_t j, const char* prefix) const;

private:
  STensorFunc* anasol; //!< Analytical stress resultant field
  VecFunc*     ana2nd; //!< Analytical 2nd derivatives of primary solution
  short int    nOrder; //!< Neumann order flag, 1 = moments, 2 = shear forces
};

#endif
