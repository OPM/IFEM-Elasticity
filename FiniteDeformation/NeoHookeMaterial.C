// $Id$
//==============================================================================
//!
//! \file NeoHookeMaterial.C
//!
//! \date Mar 08 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Neo-Hookean hyperelastic material model.
//!
//==============================================================================

#include "NeoHookeMaterial.h"
#include "FiniteElement.h"
#include "Utilities.h"
#include "Function.h"
#include "Tensor.h"
#include "IFEM.h"
#include "tinyxml2.h"


NeoHookeMaterial::NeoHookeMaterial (int ver) : LinIsotropic(false)
{
  mVER = ver / 10;
  iVOL = ver % 10;

  this->findLameParams();

  sigma_p = 0.0;
}


NeoHookeMaterial::NeoHookeMaterial (double E, double v, double d, int ver)
  : LinIsotropic(E,v,d,false)
{
  mVER = ver / 10;
  iVOL = ver % 10;

  this->findLameParams();

  sigma_p = 0.0;
}


void NeoHookeMaterial::parse (const tinyxml2::XMLElement* elem)
{
  this->LinIsotropic::parse(elem);

  if (!elem->FirstChildElement() && !Efunc && !Eaging && !nuFunc)
  {
    if (Emod >= 0.0 && utl::getAttribute(elem,"K",Emod))
      IFEM::cout <<" "<< Emod;
    if (utl::getAttribute(elem,"G",nu))
      IFEM::cout <<" "<< nu;
  }

  this->findLameParams();
}


void NeoHookeMaterial::findLameParams (double E)
{
  if (nu > 0.5)
  {
    // Assume Lame' parameters (kappa and mu are specified)
    double kappa = Bmod = E > 0.0 ? E : Emod;
    double mu    = Smod = nu;
    // Calculate Young's modulus and Poisson's ratio
    Emod = 9.0*kappa*mu/(3.0*kappa + mu);
    nu   = (1.5*kappa - mu)/(3.0*kappa + mu);
  }
  else
  {
    if (E > 0.0)
      Emod = E;
    Smod = 0.5 * Emod / (1.0 + nu);
    Bmod = nu >= 0.5 ? 0.0 : Emod / (3.0 - 6.0*nu);
    if (mVER == 1 && nu < 0.5)
      Bmod -= Smod * 2.0/3.0;
  }
}


void NeoHookeMaterial::printLog () const
{
  this->LinIsotropic::printLog();
  IFEM::cout <<"NeoHookeMaterial: mVER = "<< mVER <<" iVOL = "<< iVOL
             <<" kappa = "<< Bmod <<" mu = "<< Smod << std::endl;
}


bool NeoHookeMaterial::evaluate (Matrix& C, SymmTensor& sigma, double& U,
                                 const FiniteElement& fe, const Vec3& X,
                                 const Tensor& F, const SymmTensor& eps,
                                 char iop, const TimeDomain*,
                                 const Tensor* Fpf) const
{
  double J = F.det();
  if (J == 0.0)
  {
    std::cerr <<" *** NeoHookeMaterial::evaluate: "
              <<" Singular/zero deformation gradient\n"<< F;
    return false;
  }

  size_t ncmp = sigma.size();
  C.resize(ncmp,ncmp,true);
  sigma.leftCauchyGreen(F);

  // Evaluate the scalar stiffness function or field, if defined
  double E = -1.0;
  if (Efield)
    E = Efield->valueFE(fe);
  else if (Efunc)
    E = (*Efunc)(X);
  else if (Eaging)
    E = (*Eaging)(fe.age);

  if (E > 0.0) // Update the Lame parameters
    const_cast<NeoHookeMaterial*>(this)->findLameParams(E);

  switch (mVER) {
  case 1: // Standard hyperelastic neo-Hookean model
    U = this->stdNeoHooke(J,sigma,C);
    break;

  case 2: // Modified hyperelastic neo-Hookean model
    U = this->modNeoHooke(J,sigma,C);
    break;

  default:
    std::cerr <<" *** NeoHookeMaterial: Unknown material version ("<< mVER
              <<")."<< std::endl;
    return false;
  }

  if (iop == 2 || (iop == 3 && U == 0.0))
  {
    // Transform to 2nd Piola-Kirchhoff stresses,
    // via pull-back to reference configuration
    Tensor Fi(Fpf ? *Fpf : F);
    J = Fi.inverse();
    if (iop == 2)
    {
      sigma.transform(Fi); // sigma = F^-1 * sigma * F^-t
      sigma *= J;
      //TODO: Also pull-back the C-matrix (Total Lagrange formulation)
      std::cerr <<" *** NeoHookeMaterial::evaluate: Not available for"
                <<" Total Lagrangian formulation, sorry."<< std::endl;
      return false;
    }
    else
    {
      SymmTensor S(sigma); // sigma should be Cauchy stress when iop=3
      U = S.transform(Fi).innerProd(eps)*J;
      //TODO: Replace the above by proper path integral
      static bool first = true;
      if (first)
	std::cerr <<"  ** NeoHookeMaterial::evaluate: Path-integration of"
		  <<" strain energy density is not implemented.\n"
		  <<"     Warning: The strain energy will be incorrect."
		  << std::endl;
      first = false;
    }
  }
  else if (iop == 1) // Calculate hydrostatic pressure for output
    sigma_p = sigma.trace() / double(sigma.size() > 3 ? 3 : sigma.dim());

#ifdef INT_DEBUG
  if (iop > 0)
    std::cout <<"NeoHookeMaterial::sigma =\n"<< sigma;
  std::cout <<"NeoHookeMaterial::C ="<< C;
#endif

  return true;
}


/*!
  The energy function for the Standard Neo-Hookean material model is:
  \code
                 W(J,b) = lambda*U(J) + 0.5*mu*(I_c - 3 - 2*ln(J))
  \endcode
  where J=det(F) and U(J) is computed by the method volumetricModuli().
  The Cauchy stress tensor is:
  \code
                 Sig = [mu*(b - 1)   + lambda*Press]/J
                     =  mu*(b - 1)/J + lambda*partial_U/partial_J
  \endcode
  where b is the left Cauchy-Green deformation tensor.
*/

double NeoHookeMaterial::stdNeoHooke (double J, SymmTensor& S, Matrix& C) const
{
  double U, Press, Cv;
  this->volumetricModuli(J,U,Press,Cv);

  // Compute deviatoric and volumetric contributions
  // to the spatial stresses and constitutive tensor

  const double c1 = Smod / J;
  const double c3 = Cv * J;
  const double c2 = c1 + c1 + c3 - Press;
  const double Ic = S.trace();

  S *= c1;
  S += Press - c1;

  size_t i, j, ndim = S.dim();

  for (i = 1; i <= ndim; i++)
  {
    j = i + 3;
    C(i,i) = c2;
    C(j,j) = c1 - Press;

    for (j = 1; j <= ndim; j++)
      if (j != i)
        C(i,j) = c3 + Press;
  }

  // Compute strain energy density

  return U + Smod*(0.5*Ic - 1.5 - log(fabs(J)));
}


/*!
  Energy function for Compressible Neo-Hookean material model
  with J_2/3 regularization:
  \code
                   _ __      _            _       __
                 W(J,be) = U(J) + 0.5*mu*(J^(-2/3)be:1 - 3)
  \endcode
*/

double NeoHookeMaterial::modNeoHooke (double J, SymmTensor& S, Matrix& C) const
{
  double U, Press, Cv;
  this->volumetricModuli(J,U,Press,Cv);

  // Modify the left Cauchy-Green deformation tensor
  //                           _
  S *= pow(J,-2.0/3.0); // S = b = J^(-2/3) * b
  //                                           _
  const double I_c = S.trace() / 3.0; // = tr( b ) / 3

  size_t i, j, ndim = S.dim();

  S -= I_c;

  // Compute deviatoric part of the constitutive tensor
  //                                     _             _     _
  // Part 1: Rank one update: -2/3 G * ( b x g  +  g x b ) / J

  const double c0 = Smod * 2.0/3.0;
  const double* p = S.ptr();
  for (i = 1; i <= C.rows(); i++)
    for (j = 1; j <= ndim; j++)
    {
      C(i,j) -= c0 * p[i-1];
      C(j,i) -= c0 * p[i-1];
    }

  //                            __                     _
  // Part 2: Deviatoric term: 2 mu [ I - 1/3 g x g ] / J

  const double c1 = Smod * I_c;
  const double c2 = c1 + c1;
  const double c3 = c2 / 3.0;

  for (i = 1; i <= ndim; i++)
  {
    j = i + ndim;
    C(i,i) += c2;
    C(j,j) += c1;
    for (j = 1; j <= ndim; j++)
      C(i,j) -= c3;
  }

  // Compute deviatoric Kirchhoff stress tensor
  //           _        _
  // tau_dev = b -  tr( b ) / 3

  S *= Smod/J;
  S += Press;

  // Compute spatial deviatoric (isochoric) stresses and material moduli

  C *= 1.0/J;
  Cv *= J;

  for (i = 1; i <= ndim; i++)
  {
    C(i,i) += Cv - Press;

    j = i + ndim;
    if (j <= C.rows())
      C(j,j) -= Press;

    for (j = 1; j <= ndim; j++)
      if (j != i)
        C(i,j) += Cv + Press;
  }

  // Compute strain energy density

  return U + 1.5*Smod*(I_c - 1.0);
}


/*!
  The volumetric energy function U(J) depending on \ref iVOL as follows:

  - 1: U = K*[ 0.25*( J^2 - 1 ) - 0.5*ln(J) ]
  - 2: U = K*[ 0.5*( J - 1 )^2 ]
  - 3: U = K*[ 0.5*( ln(J) )^2 ]
  - 4: U = K*[ 2.0*( J - 1 - ln(J) ) ]

  where K = \ref Bmod.
*/

void NeoHookeMaterial::volumetricModuli (double J, double& U,
                                         double& Up, double& Upp) const
{
  switch (iVOL)
    {
    case 1: // U(J) = lambda/4 * (J^2 - 1 - 2*log(J))

      U   = 0.25 * Bmod * (J*J - 1.0 - 2.0*log(fabs(J)));
      Up  = 0.5  * Bmod * (J - 1.0/J);
      Upp = 0.5  * Bmod * (1.0 + 1.0/(J*J));
      break;

    case 2: // U(J) = lambda/2 * (J-1)^2

      U   = 0.5 * Bmod * (J-1.0)*(J-1.0);
      Up  =       Bmod * (J-1.0);
      Upp =       Bmod;
      break;

    case 3: // U(J) = lambda/2 * log(J)^2

      Upp = log(fabs(J));
      U   = 0.5 * Bmod * Upp*Upp;
      Up  =       Bmod * Upp / J;
      Upp =      (Bmod / J - Up) / J;
      break;

    case 4: // U(J) = lambda*2 * (J - 1 - log(J))

      U   = 2.0 * Bmod * (J - 1.0 - log(fabs(J)));
      Up  =       Bmod * (2.0 - 1.0/J);
      Upp =       Bmod / (J*J);
      break;

    default:
      U = Up = Upp = 0.0;
    }
}


double NeoHookeMaterial::getInternalVar (int, char* label, size_t) const
{
  if (label)
    strcpy(label,"hydrostatic pressure");

  return sigma_p;
}
