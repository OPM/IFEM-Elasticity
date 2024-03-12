// $Id$
//==============================================================================
//!
//! \file ElasticBase.C
//!
//! \date Jul 04 2014
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class representing FEM integrands for elasticity problems.
//!
//==============================================================================

#include "ElasticBase.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "TimeDomain.h"
#include "Utilities.h"
#include "BDF.h"
#include "IFEM.h"
#include "tinyxml2.h"


ElasticBase::ElasticBase () : IntegrandBase(0)
{
  nSV = 1; // Default number of solution vectors in core

  eM = eKm = eKg = 0;
  eS = gS  = iS  = dS = 0;

  memset(intPrm,0,sizeof(intPrm));

  bdf = nullptr;
}


ElasticBase::~ElasticBase ()
{
  delete bdf;
}


bool ElasticBase::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),"gravity") && nsd > 0)
  {
    utl::getAttribute(elem,"x",gravity.x);
    IFEM::cout <<"\tGravitation vector: "<< gravity.x;
    if (nsd >= 2)
    {
      utl::getAttribute(elem,"y",gravity.y);
      IFEM::cout <<" "<< gravity.y;
    }
    if (nsd >= 3)
    {
      utl::getAttribute(elem,"z",gravity.z);
      IFEM::cout <<" "<< gravity.z;
    }
    IFEM::cout << std::endl;
  }
  else
    return false;

  return true;
}


void ElasticBase::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;
  eM = eKm = eKg = 0;
  eS = gS  = iS  = 0;

  // Define element matrix/vector names, for more readable debug print.
  // Define as static, such that they don't go out of scope when exiting.

  static const char* mNames[4] = {
    "Newton", "stiffness", "mass", "geometric stiffness"
  };
  static const char* vNames[5] = {
    "external forces", "residual forces", "load gradient",
    "internal forces", "actual inertia forces"
  };

  matNames.clear();
  vecNames.clear();

  switch (mode)
    {
    case SIM::ARCLEN:
      gS  = 2;
    case SIM::STATIC:
      eKm = eKg = 1;
      eS  = iS  = 1;
      if (intPrm[3] > 0.0)
        eKg = 0; // Linear analysis, no geometric stiffness
      matNames = { mNames[0] };
      vecNames = { vNames[1], vNames[2] };
      break;

    case SIM::DYNAMIC:
      eKm = eKg = 3;
      if (intPrm[3] > 0.0)
        eKg = 0; // Linear analysis, no geometric stiffness
      else if (intPrm[3] < 0.0)
        eKg = 4; // Structural damping from material stiffness only
      eM = 2;
      eS = iS = 1;
      if (intPrm[4] == 1.0)
        eS = 3; // Store external and internal forces separately when using HHT
      matNames = { mNames[0], mNames[2], mNames[1], mNames[3] };
      vecNames = { vNames[1], vNames[4], vNames[0] };
      break;

    case SIM::BUCKLING:
      eKm = 1;
      eKg = 2;
      matNames = { mNames[1], mNames[3] };
      break;

    case SIM::VIBRATION:
      eM = 2;
    case SIM::STIFF_ONLY:
      eKm = 1;
      matNames = { mNames[1], mNames[2] };
      break;

    case SIM::MASS_ONLY:
      eM = 1;
      eS = 1;
      matNames = { mNames[2] };
      vecNames = { vNames[0] };
      break;

    case SIM::RHS_ONLY:
      eS = 1;
    case SIM::INT_FORCES:
      iS = 1;
      vecNames = { vNames[mode == SIM::RHS_ONLY ? 1 : 3] };
      break;

    default:
      ;
    }

  switch (mode)
    {
    case SIM::STATIC:
    case SIM::ARCLEN:
    case SIM::NORMS:
      primsol.resize(nSV);
      break;

    case SIM::DYNAMIC:
      primsol.resize(nSV + (intPrm[4] == 1.0 ? 4 : 2));
      break;

    case SIM::BUCKLING:
    case SIM::RHS_ONLY:
    case SIM::INT_FORCES:
    case SIM::RECOVERY:
      primsol.resize(1);
      break;

    default:
      primsol.clear();
    }
}


void ElasticBase::setIntegrationPrm (unsigned short int i, double prm)
{
  if (i < sizeof(intPrm)/sizeof(double))
    intPrm[i] = prm;
  else if (!bdf) // Using a Backward Difference Formula for time discretization
    bdf = new TimeIntegration::BDFD2(2,prm);
}


double ElasticBase::getIntegrationPrm (unsigned short int i) const
{
  return i < sizeof(intPrm)/sizeof(double) ? intPrm[i] : 0.0;
}


void ElasticBase::advanceStep (double dt, double dtn)
{
  if (bdf) bdf->advanceStep(dt,dtn);
}


std::string ElasticBase::getField1Name (size_t i, const char* prefix) const
{
  if (i > 6 || i > npv) i = 6;

  static const char* s[7] = { "u_x", "u_y", "u_z", "r_x", "r_y", "r_z",
                              "displacement" };
  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


bool ElasticBase::evalPoint (LocalIntegral& elmInt, const FiniteElement& fe,
                             const Vec3& pval)
{
  if (!eS)
  {
    std::cerr <<" *** ElasticBase::evalPoint: No load vector."<< std::endl;
    return false;
  }

  Vector& ES = static_cast<ElmMats&>(elmInt).b[eS-1];
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (unsigned short int i = 1; i <= npv && i <= 3; i++)
      ES(npv*(a-1)+i) += pval(i)*fe.N(a)*fe.detJxW;

  if (eS == 1)
  {
    RealArray& sumLoad = static_cast<ElmMats&>(elmInt).c;
    for (size_t i = 0; i < sumLoad.size() && i < 3; i++)
      sumLoad[i] += pval[i]*fe.detJxW;
  }

  return true;
}


bool ElasticBase::finalizeElement (LocalIntegral& elmInt,
                                   const TimeDomain& time, size_t)
{
  static_cast<ElmMats&>(elmInt).setStepSize(time.dt,time.it);

  return true;
}
