// $Id$
//==============================================================================
//!
//! \file KirchhoffLove.C
//!
//! \date Sep 13 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for linear Kirchhoff-Love thin plate and shell problems.
//!
//==============================================================================

#include "KirchhoffLove.h"
#include "MaterialBase.h"
#include "ElmMats.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "Tensor.h"
#include "VTF.h"


KirchhoffLove::KirchhoffLove (unsigned short int n) : IntegrandBase(n)
{
  npv = nsd < 3 ? 1 : 3; // Number of primary unknowns per node

  gravity = 0.0;
  thickness = 0.1;

  material = nullptr;
  fluxFld = nullptr;
  tracFld = nullptr;
  presFld = nullptr;
  locSys = nullptr;

  eK = eM = 0;
  eS = iS = 0;
}


KirchhoffLove::~KirchhoffLove ()
{
  if (locSys) delete locSys;
}


void KirchhoffLove::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;
  eM = eK = 0;
  eS = iS = 0;

  if (mode == SIM::RECOVERY)
    primsol.resize(1);
  else
    primsol.clear();

  switch (mode)
    {
    case SIM::STATIC:
      eK = 1;
      eS = 1;
      break;

    case SIM::VIBRATION:
      eK = 1;
      eM = 2;
      break;

    case SIM::STIFF_ONLY:
      eK = 1;
      break;

    case SIM::RHS_ONLY:
      eS = 1;
      break;

    default:
      ;
    }
}


void KirchhoffLove::initIntegration (size_t nGp, size_t nBp)
{
  presVal.clear();
  if (this->haveLoads('I'))
    presVal.resize(nGp,std::make_pair(Vec3(),Vec3()));

  tracVal.clear();
  if (this->haveLoads('B'))
    tracVal.resize(nBp,std::make_pair(Vec3(),Vec3()));
}


LocalIntegral* KirchhoffLove::getLocalIntegral (size_t nen, size_t,
                                                bool neumann) const
{
  ElmMats* result = new ElmMats();
  switch (m_mode)
    {
    case SIM::STATIC:
      result->rhsOnly = neumann;
      result->withLHS = !neumann;
      result->resize(neumann ? 0 : 1, 1);
      break;

    case SIM::VIBRATION:
      result->resize(2,0);
      break;

    case SIM::STIFF_ONLY:
      result->resize(1,0);
      break;

    case SIM::RHS_ONLY:
      result->resize(neumann ? 0 : 1, 1);

    case SIM::RECOVERY:
      result->rhsOnly = true;
      result->withLHS = false;
      break;

    default:
      ;
    }

  result->redim(npv*nen);
  return result;
}


Vec3 KirchhoffLove::getTraction (const Vec3& X, const Vec3& n) const
{
  if (fluxFld)
    return (*fluxFld)(X);
  else if (tracFld)
    return (*tracFld)(X,n);
  else
    return Vec3();
}


Vec3 KirchhoffLove::getPressure (const Vec3& X, const Vec3& n) const
{
  Vec3 p;
  p.z = material->getMassDensity(X)*gravity*thickness;

  if (presFld)
  {
    if (n.isZero())
      p.z += (*presFld)(X); // Assume pressure acts in global Z-direction
    else
      p += (*presFld)(X)*n;
  }

  return p;
}


bool KirchhoffLove::haveLoads (char type) const
{
  if (type == 'A' || type == 'I')
  {
    if (presFld)
      return true;

    if (gravity != 0.0 && material)
      if (material->getMassDensity(Vec3()) != 0.0)
        return true;
  }

  if (type == 'A' || type == 'B')
    if (fluxFld || tracFld)
      return true;

  return false;
}


void KirchhoffLove::formBodyForce (Vector& ES, const Vector& N, size_t iP,
                                   const Vec3& X, const Vec3& n,
                                   double detJW) const
{
  Vec3 p = this->getPressure(X,n);
  if (p.isZero()) return;

  if (npv == 1)
    ES.add(N,p.z*detJW);
  else for (size_t a = 1; a <= N.size(); a++)
    for (unsigned short int i = 1; i <= npv && i <= 3; i++)
      ES(npv*(a-1)+i) += N(a)*p(i)*detJW;

  // Store pressure value for visualization
  if (iP < presVal.size())
    presVal[iP] = std::make_pair(X,p);
}


void KirchhoffLove::formMassMatrix (Matrix& EM, const Vector& N,
                                    const Vec3& X, double detJW) const
{
  double rhow = material->getMassDensity(X)*thickness*detJW;
  if (rhow == 0.0) return;

  if (npv == 1)
    EM.outer_product(N,N*rhow,true);
  else
    for (size_t a = 1; a <= N.size(); a++)
      for (size_t b = 1; b <= N.size(); b++)
        for (unsigned short int i = 1; i <= npv; i++)
          EM(npv*(a-1)+i,npv*(b-1)+i) += rhow*N(a)*N(b);
}


bool KirchhoffLove::hasTractionValues () const
{
  return !tracVal.empty() || !presVal.empty();
}


bool KirchhoffLove::writeGlvT (VTF* vtf, int iStep,
                               int& geoBlk, int& nBlock) const
{
  if (tracVal.empty() && presVal.empty())
    return true;
  else if (!vtf)
    return false;

  if (!tracVal.empty())
    // Write boundary tractions as discrete point vectors to the VTF-file
    return vtf->writeVectors(tracVal,geoBlk,++nBlock,"Tractions",iStep);
  else
    // Write surface pressures as discrete point vectors to the VTF-file
    return vtf->writeVectors(presVal,geoBlk,++nBlock,"Pressure",iStep);
}


bool KirchhoffLove::evalSol (Vector& s,
                             const FiniteElement& fe, const Vec3& X,
                             const std::vector<int>& MNPC) const
{
  // Extract element displacements
  Vector eV;
  if (!primsol.empty() && !primsol.front().empty())
  {
    int ierr = utl::gather(MNPC,npv,primsol.front(),eV);
    if (ierr > 0)
    {
      std::cerr <<" *** KirchhoffLove::evalSol: Detected "
                << ierr <<" node numbers out of range."<< std::endl;
      return false;
    }
  }

  // Evaluate the stress resultants
  return this->evalSol(s,eV,fe,X,true);
}
