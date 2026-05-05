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
#include "FiniteElement.h"
#include "NewmarkMats.h"
#include "TimeDomain.h"
#include "Functions.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "Tensor.h"
#include "VTF.h"
#include "tinyxml2.h"


KirchhoffLove::KirchhoffLove (unsigned short int n, bool m) : IntegrandBase(n)
{
  npv = nsd < 3 ? 1 : 3; // Number of primary unknowns per node

  gravity = 0.0;
  constantT = true;
  thickness = new ConstFunc(0.1);

  material = nullptr;
  fluxFld = nullptr;
  tracFld = nullptr;
  linLoad = nullptr;
  locSys = nullptr;

  eK = eM = 0;
  eS = gS = iS = 0;

  includeShear = true;
  isModal = m;

  memset(intPrm,0,sizeof(intPrm));
}


KirchhoffLove::~KirchhoffLove ()
{
  if (constantT)
    delete thickness;
  delete locSys;
}


bool KirchhoffLove::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),"noshear"))
    includeShear = false;
  else if (!strcasecmp(elem->Value(),"withshear"))
    includeShear = true;
  else
    return false;

  return true;
}


void KirchhoffLove::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;
  if (mode == SIM::DYNAMIC && isModal)
    mode = m_mode = SIM::RHS_ONLY; // Modal dynamics requires RHS-vector only
  else if (mode == SIM::RHS_ONLY && !isModal && intPrm[2] != 0.0)
    mode = SIM::DYNAMIC; // RHS-only with Newmark time integration

  eM = eK = 0;
  eS = gS = iS = 0;

  if (mode == SIM::DYNAMIC)
    primsol.resize(3);
  else if (mode >= SIM::RHS_ONLY)
    primsol.resize(1);
  else
    primsol.clear();

  switch (mode)
    {
    case SIM::ARCLEN:
      gS = 2;
    case SIM::STATIC:
      eK = 1;
      eS = 1;
      break;

    case SIM::DYNAMIC:
      eK = 3;
      eM = 2;
      eS = iS = 1;
      break;

    case SIM::VIBRATION:
      eM = 2;
    case SIM::STIFF_ONLY:
      eK = 1;
      break;

    case SIM::MASS_ONLY:
      eM = 1;
      eS = 1;
      break;

    case SIM::RHS_ONLY:
      eS = 1;
    case SIM::INT_FORCES:
      iS = 1;
      break;

    default:
      ;
    }
}


/*!
  This method is overridden to optionally return the actual solution mode
  (DYNAMIC or STATIC), when the current mode flag \a m_mode is RHS_ONLY.
*/

SIM::SolutionMode KirchhoffLove::getMode (bool simMode) const
{
  if (simMode && m_mode == SIM::RHS_ONLY)
    return intPrm[2] != 0.0 ? SIM::DYNAMIC : SIM::STATIC;

  return m_mode;
}


void KirchhoffLove::setIntegrationPrm (unsigned short int i, double prm)
{
  if (i < sizeof(intPrm)/sizeof(double)) intPrm[i] = prm;
}


double KirchhoffLove::getIntegrationPrm (unsigned short int i) const
{
  return i < sizeof(intPrm)/sizeof(double) ? intPrm[i] : 0.0;
}


void KirchhoffLove::setThickness (double t)
{
  if (constantT)
    delete thickness;
  thickness = new ConstFunc(t);
  constantT = true;
}


void KirchhoffLove::setThickness (RealFunc* tf)
{
  if (constantT)
    delete thickness;
  thickness = tf;
  constantT = false;
}


void KirchhoffLove::setPressure (RealFunc* pf)
{
  if (pf)
    presFld.push_back(pf);
  else
    presFld.clear();
}


int KirchhoffLove::getIntegrandType () const
{
  int itg_type = SECOND_DERIVATIVES;
  if (m_mode == SIM::RECOVERY && includeShear)
    itg_type |= THIRD_DERIVATIVES;
  if (linLoad)
    itg_type |= INTERFACE_TERMS;

  return itg_type;
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


LocalIntegral* KirchhoffLove::getLocalIntegral (size_t nen, size_t iEl,
                                                bool neumann) const
{
  ElmMats* result = nullptr;
  if (this->inActive(iEl))
    return result; // element is not in current material group

  if (this->getMode(true) == SIM::DYNAMIC && !isModal)
    result = new NewmarkMats(intPrm[0],intPrm[1],intPrm[2],intPrm[3]);
  else
    result = new ElmMats();

  switch (m_mode)
    {
    case SIM::STATIC:
    case SIM::ARCLEN:
    case SIM::MASS_ONLY:
      result->rhsOnly = neumann;
      result->withLHS = !neumann;
      result->resize(neumann ? 0 : 1, m_mode==SIM::MASS_ONLY ? 1 : m_mode, npv);
      break;

    case SIM::DYNAMIC:
      result->rhsOnly = neumann;
      result->withLHS = !neumann;
      result->resize(neumann ? 0 : (intPrm[3] >= 0.0 ? 3 : 4),
                     neumann || intPrm[3] > 0.0 ? 1 : 2, npv);
      break;

    case SIM::VIBRATION:
      result->resize(2,0);
      break;

    case SIM::STIFF_ONLY:
      result->resize(1,0);
      break;

    case SIM::RHS_ONLY:
      if (isModal || intPrm[2] == 0.0)
        result->resize(neumann ? 0 : 1, 1, npv);
      else // RHS-only with Newmark time integration
	result->resize(neumann ? 0 : (intPrm[3] >= 0.0 ? 3 : 4),
                       neumann || intPrm[3] > 0.0 ? 1 : 2, npv);
      result->rhsOnly = true;
      result->withLHS = false;
      break;

    case SIM::RECOVERY:
      result->rhsOnly = true;
      result->withLHS = false;
      break;

    default:
      result->withLHS = false;
    }

  result->redim(npv*nen);
  return result;
}


Vec3 KirchhoffLove::getTraction (const Vec3& X, const Vec3& n, bool grd) const
{
  if (fluxFld)
    return grd ? fluxFld->timeDerivative(X) : (*fluxFld)(X);
  else if (tracFld)
    return grd ? tracFld->timeDerivative(X,n) : (*tracFld)(X,n);
  else
    return Vec3();
}


Vec3 KirchhoffLove::getPressure (const Vec3& X, const Vec3& n, bool grd) const
{
  Vec3 p;
  if (!grd)
    p.z = material->getMassDensity(X) * gravity * (*thickness)(X);

  for (RealFunc* pf : presFld)
    if (n.isZero()) // Assume pressure acts in global Z-direction
      p.z += grd ? pf->timeDerivative(X) : (*pf)(X);
    else
      p += (grd ? pf->timeDerivative(X) : (*pf)(X))*n;

  return p;
}


Vec3 KirchhoffLove::getLineLoad (const Vec3& X, const Vec3& n, bool grd) const
{
  if (!linLoad)
    return Vec3();
  else if (n.isZero()) // Assume load acts in global Z-direction
    return Vec3(0.0, 0.0, grd ? linLoad->timeDerivative(X) : (*linLoad)(X));
  else
    return (grd ? linLoad->timeDerivative(X) : (*linLoad)(X))*n;
}


bool KirchhoffLove::haveLoads (char type) const
{
  if (type == 'A' || type == 'I')
  {
    if (!presFld.empty() || linLoad)
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


void KirchhoffLove::formBodyForce (Vector& ES, RealArray& sumLoad,
                                   const Vector& N, size_t iP,
                                   const Vec3& X, const Vec3& n,
                                   double detJW, bool grd) const
{
  Vec3 p = this->getPressure(X,n,grd);
  if (p.isZero()) return;

  if (npv == 1)
    ES.add(N,p.z*detJW);
  else for (size_t a = 1; a <= N.size(); a++)
    for (unsigned short int i = 1; i <= npv; i++)
      ES(npv*(a-1)+i) += N(a)*p(i)*detJW;

  if (grd) return;

  // Integrate total external load
  if (npv == 1 && !sumLoad.empty())
    sumLoad.front() += p.z*detJW;
  else for (unsigned short int i = 0; i < npv && i < sumLoad.size(); i++)
    sumLoad[i] += p[i]*detJW;

  // Store pressure value for visualization
  if (iP < presVal.size())
    presVal[iP] = std::make_pair(X,p);
}


void KirchhoffLove::formMassMatrix (Matrix& EM, const Vector& N,
                                    const Vec3& X, double detJW) const
{
  double rhow = material->getMassDensity(X) * (*thickness)(X) * detJW;
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
  Vectors eV(primsol.size());
  for (size_t i = 0; i < primsol.size(); i++)
    if (!primsol[i].empty())
    {
      int ierr = utl::gather(MNPC,npv,primsol[i],eV[i]);
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


bool KirchhoffLove::evalPoint (LocalIntegral& elmInt, const FiniteElement& fe,
                               const Vec3& pval)
{
  if (!eS)
  {
    std::cerr <<" *** KirchhoffLove::evalPoint: No load vector."<< std::endl;
    return false;
  }

  Vector& ES = static_cast<ElmMats&>(elmInt).b[eS-1];
  if (npv == 1)
    ES.add(fe.N,pval.x*fe.detJxW);
  else for (size_t a = 1; a <= fe.N.size(); a++)
    for (unsigned short int i = 1; i <= npv; i++)
      ES(npv*(a-1)+i) += pval(i)*fe.N(a)*fe.detJxW;

  if (eS == 1)
  {
    RealArray& sumLoad = static_cast<ElmMats&>(elmInt).c;
    if (npv == 1 && !sumLoad.empty())
      sumLoad.front() += pval.x*fe.detJxW;
    else for (unsigned short int i = 0; i < npv && i < sumLoad.size(); i++)
      sumLoad[i] += pval[i]*fe.detJxW;
  }

  return true;
}


bool KirchhoffLove::finalizeElement (LocalIntegral& elmInt,
                                     const TimeDomain& time, size_t)
{
  static_cast<ElmMats&>(elmInt).setStepSize(time.dt,time.it);

  return true;
}
