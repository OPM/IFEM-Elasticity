// $Id$
//==============================================================================
//!
//! \file ElasticBar.C
//!
//! \date Aug 10 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class representing a nonlinear elastic bar.
//!
//==============================================================================

#include "ElasticBar.h"
#include "FiniteElement.h"
#include "HHTMats.h"
#include "Vec3Oper.h"
#include "IFEM.h"


ElasticBar::ElasticBar (char strm, unsigned short int nd, unsigned short int ns)
{
  strain_meassure = strm;

  npv = nd; // Number of primary unknowns per node
  nSV = ns; // Number of solution vectors in core

  stiffness = 1.0;
  lineMass = 0.0;
  lumpedMass = true;
}


void ElasticBar::printLog () const
{
  utl::LogStream& os = IFEM::cout;

  os <<"ElasticBar: Stiffness = "<< stiffness <<", Mass = "<< lineMass;
  switch (strain_meassure) {
  case 'E': os <<", Engineering strain"; break;
  case 'G': os <<", Green strain"; break;
  case 'L': os <<", Logarithmic strain"; break;
  }
  os << std::endl;
}


LocalIntegral* ElasticBar::getLocalIntegral (size_t nen, size_t, bool) const
{
  ElmMats* result;
  if (m_mode != SIM::DYNAMIC)
    result = new ElmMats();
  else if (intPrm[3] > 0.0)
    result = new NewmarkMats(intPrm[0],intPrm[1],intPrm[2],intPrm[3],
                             intPrm[4] == 2.0);
  else
    result = new HHTMats(intPrm[2],intPrm[0],intPrm[1], intPrm[4] != 1.0);

  switch (m_mode)
  {
    case SIM::STATIC:
    case SIM::MASS_ONLY:
      result->resize(1,1);
      break;

    case SIM::DYNAMIC:
      result->resize(intPrm[3] >= 0.0 ? 3 : 4,
                     intPrm[3] > 0.0 ? 1 : (intPrm[4] == 1.0 ? 3 : 2));
      break;

    case SIM::VIBRATION:
    case SIM::BUCKLING:
      result->resize(2,0);
      break;

    case SIM::STIFF_ONLY:
      result->resize(1,0);
      break;

    case SIM::RHS_ONLY:
    case SIM::INT_FORCES:
      result->resize(0,1);
      result->rhsOnly = true;
      result->withLHS = false;
      break;

    default:
      ;
  }

  result->redim(3*nen);
  return result;
}


/*!
  The integrand only depends on the position and displacement of the end points
  of the element. Thus, no basis function values or derivatives are required.
  The element stiffness matrix (including geometric terms due to the rotation
  of the bar) and mass matrix are established explicitly.
*/

bool ElasticBar::evalInt (LocalIntegral& elmInt,
                          const FiniteElement& fe,
                          const Vec3&) const
{
  // Calculate initial element length
  Vec3 X0 = fe.XC.back() - fe.XC.front();
  double L, L0 = X0.normalize();
  if (L0 <= 1.0e-8)
  {
    std::cerr <<" *** ElasticBar::evalInt: Zero initial element length "
              << L0 << std::endl;
    return false;
  }
#if INT_DEBUG > 1
  std::cout <<"ElasticBar: L0 = "<< L0;
#endif

  Vec3 X(X0), Y, Z;
  const Vector& eV = elmInt.vec.front();
  if (eV.empty())
    L = L0;
  else
  {
    // Calculate current element length
    Vec3 U1(eV.ptr(),npv);
    Vec3 U2(eV.ptr()+eV.size()-npv,npv);
    X = (fe.XC.back()+U2) - (fe.XC.front()+U1);
    L = X.normalize();
    if (L <= 1.0e-8)
    {
      std::cerr <<" *** ElasticBar::evalInt: Zero element length "
                << L << std::endl;
      return false;
    }
#if INT_DEBUG > 1
    std::cout <<" L = "<< L <<", U1 = "<< U1 <<" U2 = "<< U2;
#endif
  }

  if (eKg && L != L0)
  {
    // Determine local Y- and Z-axes
    if (fabs(X.z) > fabs(X.y))
    {
      // Define Y by projecting the global Y-axis onto the plane, then Z = X x Y
      Y.x = -X.y*X.x;
      Y.y =  X.x*X.x + X.z*X.z;
      Y.z = -X.y*X.z;
      Y.normalize();
      Z.cross(X,Y);
    }
    else
    {
      // Define Z by projecting the global Z-axis onto the plane, then Y = Z x X
      Z.x = -X.z*X.x;
      Z.y = -X.z*X.y;
      Z.z =  X.x*X.x + X.y*X.y;
      Z.normalize();
      Y.cross(Z,X);
    }
  }

  // Calculate the axial force, nominal stiffness and total mass
  double LoL0 = L/L0;
  double F = stiffness*this->getStrain(LoL0);
  double KmNom = stiffness/L0;
  double KgNom = eKg ? F/L : 0.0;
  double mass = 0.5*lineMass*L0;
  switch (strain_meassure) {
  case 'G':
    F *= LoL0;
    KmNom *= LoL0*LoL0;
    KgNom *= LoL0;
    break;
  case 'L':
    F /= LoL0;
    KmNom /= LoL0*LoL0;
    KgNom /= LoL0;
    break;
  }

#if INT_DEBUG > 1
  std::cout <<"\nElasticBar: F = "<< F <<", X = "<< X
            <<", mass = "<< mass << std::endl;
#endif

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (eKm)
  {
    Matrix& Km = elMat.A[eKm-1];
    Matrix& Kg = elMat.A[eKg ? eKg-1 : eKm-1];
    int i, j;

    // Evaluate the material stiffness matrix
    for (j = 1; j <= 3; j++)
      for (i = 1; i <= 3; i++)
        Km(i,j) = KmNom*X[i-1]*X[j-1];

    if (KgNom != 0.0)
    {
      // Evaluate the geometric stiffness matrix
      for (j = 1; j <= 3; j++)
        for (i = 1; i <= 3; i++)
        {
          double Kij = KgNom;
          switch (strain_meassure) {
          case 'E':
            Kij *= (i == j ? 1.0 : 0.0) - X[i-1]*X[j-1];
            break;
          case 'G':
            if (i != j) Kij = 0.0;
            break;
          case 'L':
            Kij *= (i == j ? 1.0 : 0.0) - 2.0*X[i-1]*X[j-1];
            break;
          default:
            Kij *= Y[i-1]*Y[j-1] + Z[i-1]*Z[j-1];
            break;
          }
          if (eKg == eKm)
            Kg(i,j) += Kij;
          else
            Kg(i,j) = Kij;
        }
    }

    for (j = 1; j <= 3; j++)
      for (i = 1; i <= 3; i++)
      {
        Km(3+i,3+j) = Km(i,j);
        Km(i,3+j) = Km(3+i,j) = -Km(i,j);
      }

    if (KgNom != 0.0 && eKg != eKm)
      for (j = 1; j <= 3; j++)
        for (i = 1; i <= 3; i++)
        {
          Kg(3+i,3+j) = Kg(i,j);
          Kg(i,3+j) = Kg(3+i,j) = -Km(i,j);
        }
  }

  if (eM)
  {
    Matrix& M = elMat.A[eM-1];
    if (lumpedMass)
      for (int i = 1; i <= 6; i++)
        M(i,i) = mass;
    else
    {
      M.fill(mass/3.0);
      for (int i = 1; i <= 6; i++)
        M(i,i) *= 2.0;
    }
  }

  if (eS && !gravity.isZero())
  {
    // External (gravitation) forces
    Vector& S = elMat.b[eS-1];
    for (int i = 1; i <= 3; i++)
      S(i) = S(3+i) = mass*gravity[i-1];
#if INT_DEBUG > 1
    std::cout <<"ElasticBar: S_ext"<< S;
#endif
  }

  if (iS)
  {
    // Internal forces
    Vector& S = elMat.b[iS-1];
    for (int i = 1; i <= 3; i++)
      if (iS == eS)
      {
        S(i)   += F*X[i-1];
        S(3+i) -= F*X[i-1];
      }
      else
      {
        S(i)   =  F*X[i-1];
        S(3+i) = -F*X[i-1];
      }
#if INT_DEBUG > 1
    if (iS == eS)
      std::cout <<"ElasticBar: S_ext - S_int"<< S;
    else
      std::cout <<"ElasticBar: -S_int"<< S;
#endif
  }

  return true;
}


double ElasticBar::getStrain (double LoverL0) const
{
  switch (strain_meassure)
  {
    case 'G': return 0.5*LoverL0*LoverL0 - 0.5;
    case 'L': return log(LoverL0);
    case 'E': return LoverL0 - 1.0;
  }

  return LoverL0 - 1.0;
}
