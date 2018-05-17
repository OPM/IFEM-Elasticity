// $Id$
//==============================================================================
//!
//! \file KirchhoffLoveShell.h
//!
//! \date Feb 25 2018
//!
//! \author Simen Skogholt Haave and Marit Gaarder Rakvaag / NTNU
//!
//! \brief Class for linear Kirchhoff-Love thin shell problems.
//!
//==============================================================================

#ifndef _KIRCHHOFF_LOVE_SHELL_H
#define _KIRCHHOFF_LOVE_SHELL_H

#include "KirchhoffLove.h"


/*!
  \brief Class representing the integrand of thin shell problems.

  \details The formulation is based on Kirchhoff-Love shell theory
  and therefore requires second-derivatives of the basis functions.
*/

class KirchhoffLoveShell : public KirchhoffLove
{
public:
  //! \brief Default constructor.
  KirchhoffLoveShell() : KirchhoffLove(3) {}
  //! \brief Empty destructor,
  virtual ~KirchhoffLoveShell() {}

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

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
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] toLocal If \e true, transform to local coordinates (if defined)
  virtual bool evalSol(Vector& s, const Vectors& eV,
                       const FiniteElement& fe, const Vec3& X,
                       bool toLocal = false) const;

  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] sm The FE in-plane stress resultant values at current point
  //! \param[out] sb The FE bending stress resultant values at current point
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] toLocal If \e true, transform to local coordinates (if defined)
  virtual bool evalSol(Vector& sm, Vector& sb, const Vectors& eV,
                       const FiniteElement& fe, const Vec3& X,
                       bool toLocal = false) const;

  //! \brief Returns max number of 2ndary solution components to print per line.
  virtual size_t getNo2ndSolPerLine() const { return 6; }

private:
  //! \brief Evaluates the stiffness matrix integrand.
  //! \param EK The element stiffness matrix to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalK(Matrix& EK, const FiniteElement& fe, const Vec3& X) const;

  //! \brief Calculates the strain-displacement matrices at current point.
  //! \param[out] Bm Membrane strain-displacement matrix
  //! \param[out] Bb Bending strain-displacement matrix
  //! \param[in] fe Finite element data at current point
  bool formBmatrix(Matrix& Bm, Matrix& Bb, const FiniteElement& fe) const;

protected:
  //! \brief Calculates some metrics at current point.
  //! \param[in] G Gradient of the cartesian coordinates at current point
  //! \param[out] g1 First covariant basis vector
  //! \param[out] g2 Second covariant basis vector
  //! \param[out] g3 Third covariant basis vector (normal vector)
  //! \param[out] n Unit normal vector
  //! \param[out] gab Covariant metric
  //! \param[out] Tlc Transformation from contravariant to local cartesian basis
  //! \return Length of the third covariant basis vector \b g3
  double getMetrics(const Matrix& G,
                    Vec3& g1, Vec3& g2, Vec3& g3, Vec3& n,
                    Vec3& gab, Matrix* Tlc = nullptr) const;

public:
  //! \brief Sets up the constitutive matrices at current point.
  //! \param[out] Dm Constitutive tensor for membrane part
  //! \param[out] Db Constitutive tensor for bending part
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] invers If \e true, the inverse matrices are establised instead
  bool formDmatrix(Matrix& Dm, Matrix& Db,
                   const FiniteElement& fe, const Vec3& X,
                   bool invers = false) const;

  //! \brief Returns the shell unit normal vector at current point.
  //! \param[in] G Gradient of the cartesian coordinates at current point
  Vec3 getShellNormal(const Matrix& G) const;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  virtual NormBase* getNormIntegrand(AnaSol*) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld = 2) const { return fld < 2 ? 4 : 18; }
  //! \brief Returns the name of the primary solution field.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField1Name(size_t i, const char* prefix) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField2Name(size_t i, const char* prefix) const;

  //! \brief Computes some derived primary solution quantities.
  virtual void primaryScalarFields(Matrix& field);
};


/*!
  \brief Class representing the integrand of Kirchhoff-Love shell energy norms.
*/

class KirchhoffLoveShellNorm : public NormBase
{
public:
  //! \brief The constructor initializes its data members.
  //! \param[in] p The thin shell problem to evaluate norms for
  KirchhoffLoveShellNorm(KirchhoffLoveShell& p);
  //! \brief Empty destructor.
  virtual ~KirchhoffLoveShellNorm() {}

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

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return true; }

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  virtual size_t getNoFields(int group) const;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  virtual std::string getName(size_t i, size_t j, const char* prefix) const;
};

#endif
