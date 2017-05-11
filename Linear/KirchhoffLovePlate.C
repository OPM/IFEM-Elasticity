// $Id$
//==============================================================================
//!
//! \file KirchhoffLovePlate.C
//!
//! \date Sep 13 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class for linear Kirchhoff-Love thin plate problems.
//!
//==============================================================================

#include "KirchhoffLovePlate.h"
#include "LinIsotropic.h"
#include "FiniteElement.h"
#include "Utilities.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "TensorFunction.h"
#include "Vec3Oper.h"
#include "AnaSol.h"
#include "VTF.h"
#include "IFEM.h"


KirchhoffLovePlate::KirchhoffLovePlate (unsigned short int n) : nsd(n)
{
  npv = 1; // Number of primary unknowns per node

  gravity = 0.0;
  thickness = 0.1;

  material = nullptr;
  locSys = nullptr;
  presFld = nullptr;
  eM = eK = 0;
  eS = 0;
}


KirchhoffLovePlate::~KirchhoffLovePlate ()
{
  if (locSys) delete locSys;
}


void KirchhoffLovePlate::printLog () const
{
  IFEM::cout <<"KirchhoffLovePlate: thickness = "<< thickness
             <<", gravity = "<< gravity << std::endl;

  if (!material)
  {
    static LinIsotropic defaultMat;
    const_cast<KirchhoffLovePlate*>(this)->material = &defaultMat;
  }

  material->printLog();
}


void KirchhoffLovePlate::setMode (SIM::SolutionMode mode)
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


LocalIntegral* KirchhoffLovePlate::getLocalIntegral (size_t nen, size_t,
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

  result->redim(nen);
  return result;
}


double KirchhoffLovePlate::getPressure (const Vec3& X) const
{
  double p = material->getMassDensity(X)*gravity*thickness;

  if (presFld)
    p += (*presFld)(X);

  return p;
}


bool KirchhoffLovePlate::haveLoads () const
{
  if (presFld) return true;

  if (gravity != 0.0 && material)
    return material->getMassDensity(Vec3()) != 0.0;

  return false;
}


void KirchhoffLovePlate::initIntegration (size_t nGp, size_t)
{
  if (this->haveLoads())
    presVal.resize(nGp,std::make_pair(Vec3(),Vec3()));
}


bool KirchhoffLovePlate::writeGlvT (VTF* vtf, int iStep,
                                    int& geoBlk, int& nBlock) const
{
  if (presVal.empty())
    return true;
  else if (!vtf)
    return false;

  // Write surface pressures as discrete point vectors to the VTF-file
  return vtf->writeVectors(presVal,geoBlk,++nBlock,"Pressure",iStep);
}


/*!
  The strain-displacement matrix for a Kirchhoff-Love plate element is formally
  defined as:
  \f[
  [B] = \left[\begin{array}{c}
   \frac{\partial^2}{\partial x^2} \\
   \frac{\partial^2}{\partial y^2} \\
  2\frac{\partial^2}{\partial x\partial y}
  \end{array}\right] [N]
  \f]
  where
  [\a N ] is the element basis functions arranged in a [NENOD] row-vector.
*/

bool KirchhoffLovePlate::formBmatrix (Matrix& Bmat,
				      const Matrix3D& d2NdX2) const
{
  const size_t nenod = d2NdX2.dim(1);
  const size_t nstrc = nsd*(nsd+1)/2;

  Bmat.resize(nstrc,nenod,true);

  if (d2NdX2.dim(2) != nsd || d2NdX2.dim(3) != nsd)
  {
    std::cerr <<" *** KirchhoffLovePlate::formBmatrix: Invalid dimension on"
	      <<" d2NdX2, "<< d2NdX2.dim(1) <<"x"<< d2NdX2.dim(2)
	      <<"x"<< d2NdX2.dim(3) <<"."<< std::endl;
    return false;
  }

  for (size_t i = 1; i <= nenod; i++)
    if (nsd == 1)
      Bmat(1,i) = d2NdX2(i,1,1);
    else
    {
      Bmat(1,i) = d2NdX2(i,1,1);
      Bmat(2,i) = d2NdX2(i,2,2);
      Bmat(3,i) = d2NdX2(i,1,2)*2.0;
    }

  return true;
}


bool KirchhoffLovePlate::formCmatrix (Matrix& C, const FiniteElement& fe,
                                      const Vec3& X, bool invers) const
{
  SymmTensor dummy(nsd); double U;
  if (!material->evaluate(C,dummy,U,fe,X,dummy,dummy, invers ? -1 : 1))
    return false;

  double factor = thickness*thickness*thickness/12.0;
  C.multiply(invers ? 1.0/factor : factor);
  return true;
}


void KirchhoffLovePlate::formMassMatrix (Matrix& EM, const Vector& N,
					 const Vec3& X, double detJW) const
{
  double rho = material->getMassDensity(X)*thickness;

  if (rho != 0.0)
    EM.outer_product(N,N*rho*detJW,true);
}


void KirchhoffLovePlate::formBodyForce (Vector& ES, const Vector& N, size_t iP,
					const Vec3& X, double detJW) const
{
  double p = this->getPressure(X);
  if (p != 0.0)
  {
    ES.add(N,p*detJW);
    // Store pressure value for visualization
    if (iP < presVal.size())
      presVal[iP] = std::make_pair(X,Vec3(0.0,0.0,p));
  }
}


bool KirchhoffLovePlate::evalInt (LocalIntegral& elmInt,
				  const FiniteElement& fe,
				  const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (eK)
  {
    // Compute the strain-displacement matrix B from d2NdX2
    Matrix Bmat;
    if (!this->formBmatrix(Bmat,fe.d2NdX2)) return false;

    // Evaluate the constitutive matrix at this point
    Matrix Cmat;
    if (!this->formCmatrix(Cmat,fe,X)) return false;

    // Integrate the stiffness matrix
    Matrix CB;
    CB.multiply(Cmat,Bmat).multiply(fe.detJxW); // CB = C*B*|J|*w
    elMat.A[eK-1].multiply(Bmat,CB,true,false,true); // EK += B^T * CB
  }

  if (eM)
    // Integrate the mass matrix
    this->formMassMatrix(elMat.A[eM-1],fe.N,X,fe.detJxW);

  if (eS)
    // Integrate the load vector due to gravitation and other body forces
    this->formBodyForce(elMat.b[eS-1],fe.N,fe.iGP,X,fe.detJxW);

  return true;
}


bool KirchhoffLovePlate::evalBou (LocalIntegral& elmInt,
				  const FiniteElement& fe,
				  const Vec3& X, const Vec3& normal) const
{
  std::cerr <<" *** KirchhoffLovePlate::evalBou not implemented."<< std::endl;
  return false;
}


bool KirchhoffLovePlate::evalSol (Vector& s,
                                  const FiniteElement& fe, const Vec3& X,
                                  const std::vector<int>& MNPC) const
{
  // Extract element displacements
  int ierr = 0;
  Vector eV;
  if (!primsol.empty() && !primsol.front().empty())
    if ((ierr = utl::gather(MNPC,1,primsol.front(),eV)))
    {
      std::cerr <<" *** KirchhoffLovePlate::evalSol: Detected "
		<< ierr <<" node numbers out of range."<< std::endl;
      return false;
    }

  // Evaluate the stress resultant tensor
  return this->evalSol(s,eV,fe,X,true);
}


bool KirchhoffLovePlate::evalSol (Vector& s, const Vector& eV,
                                  const FiniteElement& fe, const Vec3& X,
                                  bool toLocal) const
{
  if (eV.empty())
  {
    std::cerr <<" *** KirchhoffLovePlate::evalSol: No displacement vector."
	      << std::endl;
    return false;
  }
  else if (eV.size() != fe.d2NdX2.dim(1))
  {
    std::cerr <<" *** KirchhoffLovePlate::evalSol: Invalid displacement vector."
              <<"\n     size(eV) = "<< eV.size() <<"   size(d2NdX2) = "
              << fe.d2NdX2.dim(1) <<","<< fe.d2NdX2.dim(2)*fe.d2NdX2.dim(3)
              << std::endl;
    return false;
  }

  // Compute the strain-displacement matrix B from d2NdX2
  Matrix Bmat;
  if (!this->formBmatrix(Bmat,fe.d2NdX2))
    return false;

  // Evaluate the constitutive matrix at this point
  Matrix Cmat;
  if (!this->formCmatrix(Cmat,fe,X))
    return false;

  // Evaluate the curvature tensor
  SymmTensor kappa(nsd), m(nsd);
  if (!Bmat.multiply(eV,kappa)) // kappa = B*eV
    return false;

  // Evaluate the stress resultant tensor
  if (!Cmat.multiply(-1.0*kappa,m)) // m = -C*kappa
    return false;

  // Congruence transformation to local coordinate system at current point
  if (toLocal && locSys) m.transform(locSys->getTmat(X));

  s = m;
  return true;
}


size_t KirchhoffLovePlate::getNoFields (int fld) const
{
  return fld < 2 ? 1 : nsd*(nsd+1)/2;
}


std::string KirchhoffLovePlate::getField1Name (size_t,
					       const char* prefix) const
{
  if (!prefix) return "w";

  return prefix + std::string(" w");
}


std::string KirchhoffLovePlate::getField2Name (size_t i,
					       const char* prefix) const
{
  if (i >= 3) return "";

  static const char* s[3] = { "m_xx", "m_yy", "m_xy" };
  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


NormBase* KirchhoffLovePlate::getNormIntegrand (AnaSol* asol) const
{
  if (asol)
    return new KirchhoffLovePlateNorm(*const_cast<KirchhoffLovePlate*>(this),
				      asol->getStressSol());
  else
    return new KirchhoffLovePlateNorm(*const_cast<KirchhoffLovePlate*>(this));
}


KirchhoffLovePlateNorm::KirchhoffLovePlateNorm (KirchhoffLovePlate& p,
						STensorFunc* a)
  : NormBase(p), anasol(a)
{
  nrcmp = myProblem.getNoFields(2);
  projBou = true;
}


bool KirchhoffLovePlateNorm::evalInt (LocalIntegral& elmInt,
				      const FiniteElement& fe,
				      const Vec3& X) const
{
  KirchhoffLovePlate& problem = static_cast<KirchhoffLovePlate&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the inverse constitutive matrix at this point
  Matrix Cinv;
  if (!problem.formCmatrix(Cinv,fe,X,true))
    return false;

  // Evaluate the finite element stress field
  Vector mh, m, error;
  if (!problem.evalSol(mh,pnorm.vec.front(),fe,X))
    return false;

  double hk4 = fe.h*fe.h*fe.h*fe.h;
  double Res = 0.0;
  size_t ip = 0;

  // Integrate the energy norm a(w^h,w^h)
  pnorm[ip++] += mh.dot(Cinv*mh)*fe.detJxW;

  // Evaluate the body load
  double p = problem.getPressure(X);
  // Evaluate the displacement field
  double w = pnorm.vec.front().dot(fe.N);
  // Integrate the external energy (p,w^h)
  pnorm[ip++] += p*w*fe.detJxW;

#if INT_DEBUG > 3
  std::cout <<"KirchhoffLovePlateNorm::evalInt("<< fe.iel <<", "<< X <<"):";
#endif

  if (anasol)
  {
    // Evaluate the analytical stress resultant field
    m = (*anasol)(X);

    // Integrate the energy norm a(w,w)
    pnorm[ip++] += m.dot(Cinv*m)*fe.detJxW;
    // Integrate the error in energy norm a(w-w^h,w-w^h)
    error = m - mh;
    pnorm[ip++] += error.dot(Cinv*error)*fe.detJxW;

    // Residual of analytical solution (should be zero)
    if (nrcmp == 1)
      Res = anasol->dderiv(X,1,1)(1,1) + p;
    else
    {
      double d2MxxdX2 = anasol->dderiv(X,1,1)(1,1);
      double d2MyydY2 = anasol->dderiv(X,2,2)(2,2);
      double d2MxydXY = anasol->dderiv(X,1,2)(1,2);
      Res = d2MxxdX2 + d2MxydXY + d2MxydXY + d2MyydY2 + p;
    }
    // Integrate the residual error in the analytical solution
    pnorm[ip++] += hk4*Res*Res*fe.detJxW;
  }

  size_t i, j, nen = fe.N.size();
  for (i = 0; i < pnorm.psol.size(); i++)
    if (!pnorm.psol[i].empty())
    {
      // Evaluate projected stress resultant field
      Vector mr(nrcmp);
      for (j = 0; j < nrcmp; j++)
        mr[j] = pnorm.psol[i].dot(fe.N,j,nrcmp);

      // Integrate the energy norm a(w^r,w^r)
      pnorm[ip++] += mr.dot(Cinv*mr)*fe.detJxW;
      // Integrate the error in energy norm a(w^r-w^h,w^r-w^h)
      error = mr - mh;
      pnorm[ip++] += error.dot(Cinv*error)*fe.detJxW;

      // Integrate the L2-norm (m^r,m^r)
      pnorm[ip++] += mr.dot(mr)*fe.detJxW;
      // Integrate the error in L2-norm (m^r-m^h,m^r-m^h)
      pnorm[ip++] += error.dot(error)*fe.detJxW;

      // Evaluate the interior residual of the projected solution
      if (nrcmp == 1)
        Res = pnorm.psol[i].dot(fe.d2NdX2) + p;
      else
      {
        double d2mxxdx2 = pnorm.psol[i].dot(fe.d2NdX2,0,3);
        double d2myydy2 = pnorm.psol[i].dot(fe.d2NdX2,1,3,nen*3);
        double d2mxydxy = pnorm.psol[i].dot(fe.d2NdX2,2,3,nen);
        Res = d2mxxdx2 + d2mxydxy + d2mxydxy + d2myydy2 + p;
      }
      // Integrate the residual error in the projected solution
#if INT_DEBUG > 3
      std::cout <<"\n\tResidual h^4*|Laplace{m^r}-p|^2 = "
                << hk4*Res*Res << std::endl;
#endif
      pnorm[ip++] += hk4*Res*Res*fe.detJxW;

      if (anasol)
      {
        // Integrate the error in the projected solution a(w-w^r,w-w^r)
        error = m - mr;
        pnorm[ip++] += error.dot(Cinv*error)*fe.detJxW;
        ip += 2; // Make room for the local effectivity indices here
      }
    }

  return true;
}


bool KirchhoffLovePlateNorm::evalBou (LocalIntegral& elmInt,
				      const FiniteElement& fe,
				      const Vec3& X, const Vec3& normal) const
{
  std::cerr <<" *** KirchhoffLovePlateNorm::evalBou not included."<< std::endl;
  return false;
}


bool KirchhoffLovePlateNorm::finalizeElement (LocalIntegral& elmInt)
{
  if (!anasol) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate local effectivity indices as sqrt(a(e^r,e^r)/a(e,e))
  // with e^r = w^r - w^h  and  e = w - w^h,
  // and sqrt((a(e^r,e^r)+res(w^r))/a(e,e))
  for (size_t ip = 12; ip < pnorm.size(); ip += 8)
  {
    pnorm[ip-1] = sqrt(pnorm[ip-6] / pnorm[3]);
    pnorm[ip] = sqrt((pnorm[ip-6]+pnorm[ip-3]) / pnorm[3]);
  }

  return true;
}


int KirchhoffLovePlateNorm::getIntegrandType () const
{
  return SECOND_DERIVATIVES | ELEMENT_CORNERS;
}


size_t KirchhoffLovePlateNorm::getNoFields (int group) const
{
  if (group < 1)
    return this->NormBase::getNoFields();
  else if (group == 1)
    return anasol ? 5 : 2;
  else
    return anasol ? 8 : 5;
}


std::string KirchhoffLovePlateNorm::getName (size_t i, size_t j,
                                             const char* prefix) const
{
  if (i == 0 || j == 0 || j > 8 || (i == 1 && j > 5))
    return this->NormBase::getName(i,j,prefix);

  static const char* u[5] = {
    "a(w^h,w^h)^0.5",
    "(p,w^h)^0.5",
    "a(w,w)^0.5",
    "a(e,e)^0.5, e=w-w^h",
    "res(w)^0.5"
  };

  static const char* p[8] = {
    "a(w^r,w^r)^0.5",
    "a(e,e)^0.5, e=w^r-w^h",
    "(w^r,w^r)^0.5",
    "(e,e)^0.5, e=w^r-w^h",
    "res(w^r)^0.5",
    "a(e,e)^0.5, e=w-w^r",
    "effectivity index^*",
    "effectivity index^RES"
  };

  std::string name(i > 1 ? p[j-1] : u[j-1]);

  return prefix ? prefix + std::string(" ") + name : name;
}
