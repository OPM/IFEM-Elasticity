// $Id$
//==============================================================================
//!
//! \file PlasticMaterial.h
//!
//! \date Mar 16 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Elasto-plastic material models.
//!
//==============================================================================

#ifndef _PLASTIC_MATERIAL_H
#define _PLASTIC_MATERIAL_H

#include "MaterialBase.h"
#include "Tensor.h"
#include "Vec3.h"
#include <array>

class PlasticMaterial;
class ScalarFunc;


/*!
  \brief Class representing an history-dependent elasto-plastic material model.

  \details
  The plasticity models have history variables at the integration points that
  need to be "remembered" from one iteration/increment to the next. This is
  maintained inside this class as a vector of integration point objects, and
  a global counter keeping track of which integration point we are calculating.
  This counter is initialized in the beginning of each iteration, and then
  incremented by each invokation of the evaluate() method. Therefore, it is
  paramount that this method is invoked only once per integration point per
  iteration, and in the same order in every iteration. Otherwise it will
  output incorrect results.

  A separate "integration point" vector is dedicated for results points.
  This is needed because the results points (for visualization, etc.) are not
  the same as the integration points used in the tangent evaluation.
*/

class PlasticMaterial : public Material
{
  //! \brief Class representing an elasto-plastic material point.
  class PlasticPoint
  {
  public:
    //! \brief Constructor initializing the material parameters.
    //! \param[in] prm Pointer to actual material model object.
    //! \param[in] n Number of space dimensions
    PlasticPoint(const PlasticMaterial* prm, unsigned short int n);
    //! \brief Empty destructor.
    virtual ~PlasticPoint() {}

    //! \brief Updates the internal history variables after convergence.
    bool updateState(bool updateVars = true);

    //! \brief Evaluates the constitutive relation at this point.
    //! \param[out] C Constitutive matrix at current point
    //! \param[out] sigma Stress tensor at current point
    //! \param[in] F Deformation gradient at current point
    //! \param[in] prm Nonlinear solution algorithm parameters
    int evaluate(Matrix& C, SymmTensor& sigma,
                 const Tensor& F, const TimeDomain& prm) const;

    //! \brief Updates the path integral of the strain energy density.
    //! \param[in] S Stress tensor at current configuration
    //! \param[in] E Strain tensor at current configuration
    //! \return Updated strain energy density
    double energyIntegral(const SymmTensor& S, const SymmTensor& E);

    //! \brief Returns a history variable.
    double getVariable(size_t i) const { return HVc[i]; }

    //! \brief Returns whether this material point has diverged or not.
    virtual bool diverged() const { return updated == 'd'; }

  protected:
    //! \brief Evaluates the yield function and its derivatives.
    bool yfunc(bool lIter, double Epp, double I1, double J2, double J3,
               double& f1, double& f2, double& f3, double& f11, double& f22,
               double& f33, double& f12, double& f13, double& f23,
               double& Ypr, double& YY, double& Yield) const;

  private:
    const RealArray& pMAT; //!< Material property parameters
    const ScalarFunc* hfn; //!< Isotropic hardening function

    using HistoryVars = std::array<double,10>; //!< History variable container

    HistoryVars HVc; //!< History variable values in current configuration
    HistoryVars HVp; //!< History variable values in previous configuration

    char updated; //!< Flag indicating whether history variables are updated

    // Data for path integral of strain energy
    SymmTensor Ep; //!< Strain tensor, previous configuration
    SymmTensor Sp; //!< Stress tensor, previous configuration
    double     Up; //!< Strain energy density

  public:
    Tensor Fp; //!< Deformation gradient, previous configuration
  };

  //! \brief Class representing an elasto-plastic result point.
  class ResultPoint : public PlasticPoint
  {
  public:
    //! \brief Constructor initializing the material parameters.
    //! \param[in] prm Pointer to actual material model object.
    //! \param[in] n Number of space dimensions
    ResultPoint(const PlasticMaterial* prm, unsigned short int n);

    //! \brief Evaluates principal stresses at this point.
    //! \param[in] sigma Stress tensor at current point
    void principalStress(const SymmTensor& sigma);

    //! \brief Returns the mean stress:
    double getMeanStress() const { return sigma_h; }
    //! \brief Returns the principal stress.
    double getPrnStress(size_t i = 0) const { return prin[i]; }
    //! \brief Returns the stress triaxiality.
    double getTriax() const;
    //! \brief Returns the Lode parameter.
    double getLode() const;

  private:
    double sigma_h; //!< Mean stress
    double sigma_e; //!< von Mises stress
    Vec3   prin;    //!< Principal stresses
  };

public:
  //! \brief Constructor initializing the material parameters.
  PlasticMaterial(const RealArray& p, const ScalarFunc* hcurve = nullptr);
  //! \brief The destructor frees the dynamically allocated data.
  virtual ~PlasticMaterial();

  //! \brief Prints out material parameters to the log stream.
  virtual void printLog() const;

  //! \brief Initializes the material with the number of integration points.
  virtual void initIntegration(size_t nGP);
  //! \brief Initializes the material model for a new integration loop.
  //! \param[in] prm Nonlinear solution algorithm parameters
  //!
  //! \details This method is invoked once before starting the numerical
  //! integration over the entire spatial domain, and is used to reset the
  //! global integration point counter, and to update the history variables.
  virtual void initIntegration(const TimeDomain& prm);
  //! \brief Initializes the material model for a new result point loop.
  virtual void initResultPoints();

  //! \brief Evaluates the mass density at current point.
  virtual double getMassDensity(const Vec3&) const { return pMAT[3]; }

  using Material::evaluate;
  //! \brief Evaluates the constitutive relation at current integration point.
  //! \param[out] C Constitutive matrix at current point
  //! \param[out] sigma Stress tensor at current point
  //! \param[out] U Strain energy density at current point
  //! \param[in] fe Finite element quantities at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] F Deformation gradient at current point
  //! \param[in] eps Strain tensor at current point
  //! \param[in] iop Calculation option:
  //! - 0 : Calculate the constitutive matrix only
  //! - 1 : Cauchy stresses and the tangent constitutive matrix
  //! - 2 : 2nd Piola-Kirchhoff stresses and the tangent constitutive matrix
  //! - 3 : Calculate strain energy density only
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] Fpf Deformation gradient for push-forward transformation
  virtual bool evaluate(Matrix& C, SymmTensor& sigma, double& U,
                        const FiniteElement& fe, const Vec3& X,
                        const Tensor& F, const SymmTensor& eps,
                        char iop, const TimeDomain* prm,
                        const Tensor* Fpf = nullptr) const;

  //! \brief Returns whether the material model has diverged.
  //! \param[in] iP1 Global (1-based) index for the integration point to check,
  //! checking all points if \a iP1 is zero
  virtual bool diverged(size_t iP1) const;

  //! \brief Returns number of internal result variables of the material model.
  virtual int getNoIntVariables() const { return 7; }

  //! \brief Returns an internal variable associated with the material model.
  //! \param[in] idx Index (1-based) of the internal variable
  //! \param[out] label Name of the internal variable (for result presentation)
  //! \param[in] iP1 Global (1-based) index for current point
  virtual double getInternalVar(int idx, char* label, size_t iP1) const;

private:
  RealArray pMAT; //!< Material property parameters

  const ScalarFunc* hardening; //!< Isotropic hardening function

  bool iAmIntegrating; //!< Flag indicating integration or result evaluation
  bool firstStep;      //!< \e true when we are doing the first time step

  mutable size_t                    iP2;       //!< Global result point counter
  mutable std::vector<ResultPoint*> itgPoints; //!< Integration point data
  mutable std::vector<ResultPoint*> resPoints; //!< Result point data
};

#endif
