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
#include "Tensor.h"
#include "VTF.h"


KirchhoffLove::KirchhoffLove (unsigned short int n) : IntegrandBase(n)
{
  npv = nsd < 3 ? 1 : 3; // Number of primary unknowns per node

  gravity = 0.0;
  thickness = 0.1;

  material = nullptr;
  locSys = nullptr;
  presFld = nullptr;
  eM = eK = 0;
  eS = 0;
}


KirchhoffLove::~KirchhoffLove ()
{
  if (locSys) delete locSys;
}


void KirchhoffLove::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;
  eM = eK = 0;
  eS = 0;

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


double KirchhoffLove::getPressure (const Vec3& X) const
{
  double p = material->getMassDensity(X)*gravity*thickness;

  if (presFld)
    p += (*presFld)(X);

  return p;
}


bool KirchhoffLove::haveLoads () const
{
  if (presFld)
    return true;

  if (gravity != 0.0 && material)
    return material->getMassDensity(Vec3()) != 0.0;

  return false;
}


void KirchhoffLove::initIntegration (size_t nGp, size_t)
{
  if (this->haveLoads())
    presVal.resize(nGp,std::make_pair(Vec3(),Vec3()));
}


bool KirchhoffLove::writeGlvT (VTF* vtf, int iStep,
                               int& geoBlk, int& nBlock) const
{
  if (presVal.empty())
    return true;
  else if (!vtf)
    return false;

  // Write surface pressures as discrete point vectors to the VTF-file
  return vtf->writeVectors(presVal,geoBlk,++nBlock,"Pressure",iStep);
}
