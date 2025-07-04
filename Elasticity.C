// $Id$
//==============================================================================
//!
//! \file Elasticity.C
//!
//! \date Nov 12 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for linear and nonlinear elasticity problems.
//!
//==============================================================================

#include "Elasticity.h"
#include "ElasticityUtils.h"
#include "GlobalIntegral.h"
#include "IsotropicTextureMat.h"
#include "FiniteElement.h"
#include "HHTMats.h"
#include "BDFMats.h"
#include "ElmNorm.h"
#include "AnaSol.h"
#include "TensorFunction.h"
#include "Vec3Oper.h"
#include "Utilities.h"
#include "VTF.h"
#include "IFEM.h"
#include "tinyxml2.h"
#include <iomanip>

#ifndef epsR
//! \brief Zero tolerance for the radial coordinate.
#define epsR 1.0e-16
#endif

bool Elasticity::wantStrain          = false;
bool Elasticity::wantPrincipalStress = false;
bool Elasticity::asolProject         = false;


Elasticity::Elasticity (unsigned short int n, bool ax) : dS(0), axiSymmetry(ax)
{
  nsd = axiSymmetry ? 2 : n;
  nDF = axiSymmetry ? 3 : nsd;
  npv = nsd; // Number of primary unknowns per node

  // Assign default material properties, in case of no user-input
  static LinIsotropic defaultMat;
  material = &defaultMat;

  locSys  = nullptr;
  tracFld = nullptr;
  fluxFld = nullptr;
  bodyFld = nullptr;
  pDirBuf = nullptr;
  dualRHS = nullptr;
  myReacI = nullptr;

  gamma = 1.0;

  calcMaxVal = true;
}


Elasticity::~Elasticity ()
{
  if (locSys) delete locSys;
  if (pDirBuf) delete pDirBuf;
}


bool Elasticity::parse (const tinyxml2::XMLElement* elem)
{
  if (this->ElasticBase::parse(elem))
    return true;
  else if (!strcasecmp(elem->Value(),"stabilization"))
  {
    utl::getAttribute(elem,"gamma",gamma);
    IFEM::cout <<"\tStabilization parameter "<< gamma << std::endl;
  }
  else if (!strcasecmp(elem->Value(),"localsystem"))
    this->parseLocalSystem(elem);
  else
    return false;

  return true;
}


Material* Elasticity::parseMatProp (char* cline)
{
  bool planeStress = nsd == 2 ? !Elastic::planeStrain : false;

  double E   = atof(strtok(cline," "));
  double nu  = atof(strtok(nullptr," "));
  double rho = atof(strtok(nullptr," "));
  IFEM::cout << E <<" "<< nu <<" "<< rho << std::endl;
  material = new LinIsotropic(E,nu,rho,planeStress,axiSymmetry);
  return material;
}


Material* Elasticity::parseMatProp (const tinyxml2::XMLElement* elem)
{
  utl::getAttribute(elem,"planeStrain",Elastic::planeStrain);
  bool planeStress = nsd == 2 ? !Elastic::planeStrain : false;

  if (!strcasecmp(elem->Value(),"texturematerial"))
    material = new IsotropicTextureMat(planeStress,axiSymmetry);
  else
    material = new LinIsotropic(planeStress,axiSymmetry);
  material->parse(elem);
  return material;
}


void Elasticity::printLog () const
{
  utl::LogStream& os = IFEM::cout;

  if (axiSymmetry)
    os <<"Axial-symmetric Elasticity problem\n";
  os <<"Elasticity: "<< nsd <<"D, gravity =";
  for (unsigned short int d = 0; d < nsd; d++)
    os <<" "<< gravity[d];
  os << std::endl;

  material->printLog();
}


void Elasticity::setMaterial (Material* mat)
{
  if (mat == material) return;

#ifdef INT_DEBUG
  if (mat && material)
  {
    IFEM::cout <<"\nElasticity::setMaterial: Switching material properties:\n";
    mat->printLog();
  }
#endif
  material = mat;
}


void Elasticity::addExtrFunction (FunctionBase* extr)
{
  if (dynamic_cast<VecFunc*>(extr))
  {
    dualFld.push_back(extr);
    if (!dualRHS) dualRHS = extr;
    // Set number of solution states
    nSV = nCS = 1+dualFld.size();
    dS = 1;
  }
  else
  {
    dualFld.clear();
    dualRHS = nullptr;
  }
}


LocalIntegral* Elasticity::getLocalIntegral (size_t nen, size_t iEl,
					     bool neumann) const
{
  ElmMats* result = nullptr;
  if (this->inActive(iEl))
    return result; // element is not in current material group

  const bool linDyn = intPrm[3] > 0.0;
  const bool useHHT = intPrm[4] == 1.0;
  const bool useGA  = intPrm[4] == 2.0;

  if (bdf) // time integration using Backward Difference Formula
    result = new BDFMats(*bdf);
  else if (this->getMode(true) != SIM::DYNAMIC)
    result = new ElmMats(); // linear or nonlinear (quasi-)static analysis
  else if (linDyn) // linear dynamic analysis with Newmark time integration
    result = new NewmarkMats(intPrm[0], intPrm[1], intPrm[2], intPrm[3], useGA);
  else // nonlinear dynamic analysis with Newmark time integration
    result = new HHTMats(intPrm[2], intPrm[0], intPrm[1], !useHHT);

  switch (bdf && m_mode == SIM::STATIC ? SIM::DYNAMIC : m_mode)
  {
    case SIM::STATIC:
    case SIM::MASS_ONLY:
      result->rhsOnly = neumann;
      result->withLHS = !neumann;
      result->resize(neumann ? 0 : 1, neumann || dS == 0 ? 1 : dS, nsd);
      break;

    case SIM::ARCLEN:
      result->rhsOnly = neumann;
      result->withLHS = !neumann;
      result->resize(neumann ? 0 : 1, 2, nsd);
      break;

    case SIM::DYNAMIC:
      result->rhsOnly = neumann;
      result->withLHS = !neumann;
      result->resize(neumann ? 0 : (intPrm[3] < 0.0 ? 4 : 3),
                     useHHT ? 3 : (neumann || linDyn ? 1 : 2), nsd);
      break;

    case SIM::VIBRATION:
    case SIM::BUCKLING:
      result->resize(2,0);
      break;

    case SIM::STIFF_ONLY:
      result->resize(1,0);
      break;

    case SIM::RHS_ONLY:
      if (nSV > nCS) // RHS-only with Newmark time integration
        result->resize(neumann ? 0 : (intPrm[3] < 0.0 ? 4 : 3),
                       useHHT ? 3 : (neumann || linDyn ? 1 : 2), nsd);
      else
        result->resize(neumann ? 0 : 1, 1, nsd);
      result->rhsOnly = true;
      result->withLHS = false;
      break;

    case SIM::INT_FORCES:
      result->resize(neumann ? 0 : 1, 1, nsd);

    case SIM::RECOVERY:
      result->rhsOnly = true;
      result->withLHS = false;
      break;

    default:
      result->withLHS = false;
  }

  result->redim(nsd*nen);
  result->Aname = matNames;
  result->Bname = vecNames;
  return result;
}


bool Elasticity::hasBoundaryTerms () const
{
  return (m_mode >= SIM::STATIC && m_mode <= SIM::DYNAMIC) ||
    (m_mode >= SIM::RHS_ONLY && m_mode <= SIM::INT_FORCES) ||
    m_mode == SIM::NORMS;
}


void Elasticity::setSecondaryInt (GlobalIntegral* gq)
{
  delete myReacI;
  myReacI = gq;
}


GlobalIntegral& Elasticity::getGlobalInt (GlobalIntegral* gq) const
{
  if (m_mode == SIM::RHS_ONLY && myReacI)
    return *myReacI;

  return this->ElasticBase::getGlobalInt(gq);
}


Vec3 Elasticity::getTraction (const Vec3& X, const Vec3& n, bool grd) const
{
  if (fluxFld)
    return grd ? fluxFld->timeDerivative(X) : (*fluxFld)(X);
  else if (tracFld)
    return grd ? tracFld->timeDerivative(X,n) : (*tracFld)(X,n);
  else
    return Vec3();
}


Vec3 Elasticity::getBodyforce (const Vec3& X, bool grd) const
{
  if (grd)
    return bodyFld ? bodyFld->deriv(X,4) : Vec3();

  Vec3 f(gravity);
  f *= material->getMassDensity(X);

  if (bodyFld)
    f += (*bodyFld)(X);

  return f;
}


bool Elasticity::haveLoads () const
{
  if (tracFld) return true;
  if (fluxFld) return true;
  if (bodyFld) return true;

  for (unsigned short int i = 0; i < nsd; i++)
    if (gravity[i] != 0.0)
      return material->getMassDensity(Vec3()) != 0.0;

  return false;
}


void Elasticity::initIntegration (size_t, size_t nBp)
{
  tracVal.clear();
  tracVal.resize(nBp,std::make_pair(Vec3(),Vec3()));
}


void Elasticity::initResultPoints (double, bool prinDir)
{
  if (wantPrincipalStress && prinDir)
  {
    if (!pDirBuf) pDirBuf = new Vec3Vec();
    pDirBuf->clear();
  }
  else if (pDirBuf)
  {
    delete pDirBuf;
    pDirBuf = nullptr;
  }
}


bool Elasticity::writeGlvT (VTF* vtf, int iStep, int& geoBlk, int& nBlock) const
{
  if (tracVal.empty())
    return true;
  else if (!vtf)
    return false;

  // Write boundary tractions as discrete point vectors to the VTF-file
  return vtf->writeVectors(tracVal,geoBlk,++nBlock,"Tractions",iStep);
}


/*!
  The strain-displacement matrix for a continuum element is formally defined as:
  \f[ \mbox{In 3D,~~}
  [B] = \left[\begin{array}{ccc}
  \frac{\partial}{\partial x} & 0 & 0 \\
  0 & \frac{\partial}{\partial y} & 0 \\
  0 & 0 & \frac{\partial}{\partial z} \\
  \frac{\partial}{\partial y} & \frac{\partial}{\partial x} & 0 \\
  0 & \frac{\partial}{\partial z} & \frac{\partial}{\partial y} \\
  \frac{\partial}{\partial z} & 0 & \frac{\partial}{\partial x}
  \end{array}\right] [N] \hskip 5mm\mbox{and in 2D,~~}
  [B] = \left[\begin{array}{ccc}
  \frac{\partial}{\partial x} & 0 \\
  0 & \frac{\partial}{\partial y} \\
  \frac{\partial}{\partial y} & \frac{\partial}{\partial x}
  \end{array}\right] [N]
  \f]
  where
  [\a N ] is the element basis functions arranged in a [nsd][nsd*NENOD] matrix.
*/

bool Elasticity::formBmatrix (Matrix& Bmat, const Matrix& dNdX) const
{
  const size_t nenod = dNdX.rows();
  const size_t nstrc = nsd*(nsd+1)/2;
  Bmat.resize(nstrc*nsd,nenod,true);
  if (dNdX.cols() < nsd)
  {
    std::cerr <<" *** Elasticity::formBmatrix: Invalid dimension on dNdX, "
	      << dNdX.rows() <<"x"<< dNdX.cols() <<"."<< std::endl;
    return false;
  }

#define INDEX(i,j) i+nstrc*(j-1)

  switch (nsd) {
  case 1:

    // Strain-displacement matrix for 1D elements:
    //
    //   [B] = | d/dx | * [N]

    for (size_t i = 1; i <= nenod; i++)
      Bmat(1,i) = dNdX(i,1);
    break;

  case 2:

    // Strain-displacement matrix for 2D elements:
    //
    //         | d/dx   0   |
    //   [B] = |  0    d/dy | * [N]
    //         | d/dy  d/dx |

    for (size_t i = 1; i <= nenod; i++)
    {
      // Normal strain part
      Bmat(INDEX(1,1),i) = dNdX(i,1);
      Bmat(INDEX(2,2),i) = dNdX(i,2);
      // Shear strain part
      Bmat(INDEX(3,1),i) = dNdX(i,2);
      Bmat(INDEX(3,2),i) = dNdX(i,1);
    }
    break;

  case 3:

    // Strain-displacement matrix for 3D elements:
    //
    //         | d/dx   0     0   |
    //         |  0    d/dy   0   |
    //   [B] = |  0     0    d/dz | * [N]
    //         | d/dy  d/dx   0   |
    //         |  0    d/dz  d/dy |
    //         | d/dz   0    d/dx |

    for (size_t i = 1; i <= nenod; i++)
    {
      // Normal strain part
      Bmat(INDEX(1,1),i) = dNdX(i,1);
      Bmat(INDEX(2,2),i) = dNdX(i,2);
      Bmat(INDEX(3,3),i) = dNdX(i,3);
      // Shear strain part
      Bmat(INDEX(4,1),i) = dNdX(i,2);
      Bmat(INDEX(4,2),i) = dNdX(i,1);
      Bmat(INDEX(5,2),i) = dNdX(i,3);
      Bmat(INDEX(5,3),i) = dNdX(i,2);
      Bmat(INDEX(6,1),i) = dNdX(i,3);
      Bmat(INDEX(6,3),i) = dNdX(i,1);
    }
    break;

  default:
    std::cerr <<" *** Elasticity::formBmatrix: nsd="<< nsd << std::endl;
    return false;
  }

#undef INDEX

  Bmat.resize(nstrc,nsd*nenod);
  return true;
}


/*!
  The strain-displacement matrix for an axially symmetric 3D continuum element
  is formally defined as:
  \f[
  [B] = \left[\begin{array}{cc}
  \frac{\partial}{\partial r} &                0            \\
                 0            & \frac{\partial}{\partial z} \\
         \frac{1}{r}          &                0            \\
  \frac{\partial}{\partial z} & \frac{\partial}{\partial r}
  \end{array}\right] [N]
  \f]
  where
  [\a N ] is the element basis functions arranged in a [2][2*NENOD] matrix.
*/

bool Elasticity::formBmatrix (Matrix& Bmat, const Vector& N, const Matrix& dNdX,
			      const double r) const
{
  const size_t nenod = N.size();
  Bmat.resize(8,nenod,true);
  if (dNdX.cols() < 2)
  {
    std::cerr <<" *** Elasticity::formBmatrix: Invalid dimension on dNdX, "
	      << dNdX.rows() <<"x"<< dNdX.cols() <<"."<< std::endl;
    return false;
  }
  else if (r < -epsR)
  {
    std::cerr <<" *** Elasticity::formBmatrix: Invalid point r < 0, "
	      << r << std::endl;
    return false;
  }

#define INDEX(i,j) i+4*(j-1)

  // Strain-displacement matrix for 3D axisymmetric elements:
  //
  //         | d/dr   0   |
  //   [B] = |  0    d/dz | * [N]
  //         | 1/r    0   |
  //         | d/dz  d/dr |

  for (size_t i = 1; i <= nenod; i++)
  {
    // Normal strain part
    Bmat(INDEX(1,1),i) = dNdX(i,1);
    Bmat(INDEX(2,2),i) = dNdX(i,2);
    // Hoop strain part
    Bmat(INDEX(3,1),i) = r <= epsR ? dNdX(i,1) : N(i)/r;
    // Shear strain part
    Bmat(INDEX(4,1),i) = dNdX(i,2);
    Bmat(INDEX(4,2),i) = dNdX(i,1);
  }

#undef INDEX

  Bmat.resize(4,2*nenod);
  return true;
}


bool Elasticity::formDefGradient (const Vector& eV, const Vector& N,
                                  const Matrix& dNdX, double r, Tensor& F,
                                  bool gradOnly) const
{
  F = gradOnly ? 0.0 : 1.0;
  if (eV.empty())
    return true; // Initial state, deformation gradient is identity tensor

  const size_t nenod = dNdX.rows();
  if (eV.size() != nsd*nenod || dNdX.cols() < nsd)
  {
    std::cerr <<" *** Elasticity::formDefGradient: Invalid dimension,"
              <<" dNdX("<< nenod <<","<< dNdX.cols() <<")"<< std::endl;
    return false;
  }

  // Compute the deformation gradient, [F] = [I] + [dudX] = [I] + [dNdX]*[u].
  // Notice that the matrix multiplication method used here treats the element
  // displacement vector, *eV, as a matrix whose number of columns equals the
  // number of rows in the matrix dNdX.
  Matrix dUdX;
  if (!dUdX.multiplyMat(eV,dNdX)) // dUdX = Grad{u} = eV*dNdX
    return false;

  // Cannot use operator= here, in case F is of higher dimension than dUdX
  for (size_t i = 1; i <= dUdX.rows(); i++)
    for (size_t j = 1; j <= dUdX.cols(); j++)
      F(i,j) += dUdX(i,j);

  // Add the dU/r term to the F(3,3)-term for axisymmetric problems
  if (axiSymmetry && r > epsR && !gradOnly)
    F(3,3) += eV.dot(N,0,nsd)/r;

#if INT_DEBUG > 0
  std::cout <<"Elasticity::eV ="<< eV
            <<"Elasticity::F =\n"<< F;
#endif

  return true;
}


bool Elasticity::kinematics (const Vector& eV,
			     const Vector& N, const Matrix& dNdX, double r,
			     Matrix& B, Tensor&, SymmTensor& eps) const
{
  // Evaluate the strain-displacement matrix, B
  if (axiSymmetry)
  {
    if (!this->formBmatrix(B,N,dNdX,r))
      return false;
  }
  else
  {
    if (!this->formBmatrix(B,dNdX))
      return false;
  }

  if (eV.empty() || eps.dim() == 0)
    return true;

  // Evaluate the strains
  return B.multiply(eV,eps); // eps = B*eV
}


void Elasticity::formKG (Matrix& EM, const Vector& N, const Matrix& dNdX,
			 double r, const Tensor& sigma, double detJW) const
{
#if INT_DEBUG > 3
  std::cout <<"Elasticity::sigma =\n"<< sigma
            <<"Elasticity::kg =";
#endif

  unsigned short int i, j;
  double kgrr = axiSymmetry && r > 0.0 ? sigma(3,3)/(r*r) : 0.0;
  for (size_t a = 1; a <= dNdX.rows(); a++)
    for (size_t b = 1; b <= dNdX.rows(); b++)
    {
      double kg = 0.0;
      for (i = 1; i <= nsd; i++)
	for (j = 1; j <= nsd; j++)
	  kg += dNdX(a,i)*sigma(i,j)*dNdX(b,j);
#if INT_DEBUG > 3
      std::cout << (b == 1 ? '\n' : ' ') << kg;
#endif

      for (i = 1; i <= nsd; i++)
	EM(nsd*(a-1)+i,nsd*(b-1)+i) += kg*detJW;

      if (kgrr > 0.0)
	EM(nsd*(a-1)+1,nsd*(b-1)+1) += N(a)*kgrr*N(b)*detJW;
    }
#if INT_DEBUG > 3
  std::cout << std::endl;
#endif
}


void Elasticity::formMassMatrix (Matrix& EM, const Vector& N,
				 const Vec3& X, double detJW) const
{
  double rhow = material->getMassDensity(X)*detJW;
  if (rhow == 0.0) return;

  for (size_t a = 1; a <= N.size(); a++)
    for (size_t b = 1; b <= N.size(); b++)
      for (unsigned short int i = 1; i <= nsd; i++)
	EM(nsd*(a-1)+i,nsd*(b-1)+i) += rhow*N(a)*N(b);
}


void Elasticity::formBodyForce (Vector& ES, RealArray& sumLoad, const Vector& N,
                                const Vec3& X, double detJW, bool grd) const
{
  Vec3 f = this->getBodyforce(X,grd);
  if (f.isZero()) return;

  f *= detJW;
  for (size_t a = 1; a <= N.size(); a++)
    for (unsigned short int i = 1; i <= nsd; i++)
      ES(nsd*(a-1)+i) += f[i-1]*N(a);

  if (grd) return;

  // Integrate total external load
  for (unsigned short int i = 0; i < nsd && i < sumLoad.size(); i++)
    sumLoad[i] += f[i];
}


bool Elasticity::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                          const Vec3& X, const Vec3& normal) const
{
  if (!tracFld && !fluxFld)
  {
    std::cerr <<" *** Elasticity::evalBou: No tractions."<< std::endl;
    return false;
  }
  else if (!eS)
  {
    std::cerr <<" *** Elasticity::evalBou: No load vector."<< std::endl;
    return false;
  }

  // Axi-symmetric integration point volume; 2*pi*r*|J|*w
  const double detJW = axiSymmetry ? 2.0*M_PI*X.x*fe.detJxW : fe.detJxW;

  // Evaluate the surface traction
  Vec3 T = this->getTraction(X,normal);

  // Store traction value for visualization
  if (fe.iGP < tracVal.size())
  {
    tracVal[fe.iGP].first = X;
    tracVal[fe.iGP].second += T;
  }

  // Pull-back traction to reference configuration
  if (!this->pullBackTraction(T))
    return false;

  // Integrate the force vector
  Vector& ES = static_cast<ElmMats&>(elmInt).b[eS-1];
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (unsigned short int i = 1; i <= nsd; i++)
      ES(nsd*(a-1)+i) += T[i-1]*fe.N(a)*detJW;

  // Integrate total external load
  RealArray& sumLoad = static_cast<ElmMats&>(elmInt).c;
  for (unsigned short int i = 0; i < nsd && i < sumLoad.size(); i++)
    sumLoad[i] += T[i]*fe.detJxW;

  if (gS)
  {
    // Evaluate the surface traction gradient
    T = this->getTraction(X,normal,true);

    // Pull-back traction to reference configuration
    if (!this->pullBackTraction(T))
      return false;

    // Integrate the force gradient vector
    Vector& GS = static_cast<ElmMats&>(elmInt).b[gS-1];
    for (size_t a = 1; a <= fe.N.size(); a++)
      for (unsigned short int i = 1; i <= nsd; i++)
        GS(nsd*(a-1)+i) += T[i-1]*fe.N(a)*detJW;
  }

  return true;
}


Vec3 Elasticity::evalSol (const Vector& eV, const Vector& N) const
{
  Vec3 u;
  if (eV.size() == nsd*N.size())
    for (unsigned short int i = 0; i < nsd; i++)
      u[i] = eV.dot(N,i,nsd);

  return u;
}


bool Elasticity::formCmat (Matrix& C, const FiniteElement& fe,
                           const Vec3& X, bool inverted) const
{
  SymmTensor dummy(nsd,axiSymmetry); double U;
  return material->evaluate(C,dummy,U,fe,X,dummy,dummy, inverted ? -1 : 0);
}


bool Elasticity::evalSol2 (Vector& s, const Vectors& eV,
                           const FiniteElement& fe, const Vec3& X) const
{
  Vec3* pBuf = nullptr;
  if (pDirBuf)
  {
    // Store principal stress directions in the internal buffer
    size_t ifirst = pDirBuf->size();
    pDirBuf->resize(ifirst+2);
    pBuf = pDirBuf->data() + ifirst;
  }

  // Evaluate the stress tensor
  if (fe.detJxW == 0.0)
  {
    // Singular point, just return an empty vector for now
    s.clear();
    return true;
  }
  else if (!this->evalSol(s,eV,fe,X,true,pBuf))
    return false;
#if INT_DEBUG > 2
  else if (pBuf)
    std::cout <<"Elasticity::evalSol2("<< X <<"): "
	      <<" Pdir1 = "<< pBuf[0] <<", Pdir2 = "<< pBuf[1] << std::endl;
#endif

  // Additional result variables?
  for (int i = 1; i <= material->getNoIntVariables(); i++)
    s.push_back(material->getInternalVariable(i,nullptr,fe.iGP));

  if (!calcMaxVal || maxVal.empty())
    return true; // Avoid thread sync if no max value calculation

  // Find the maximum values for each quantity. This block must be performed
  // serially on multi-threaded runs too, due to the update of the maxVal array
  // which is a member of the Elasticity class. Therefore the critical pragma.
#pragma omp critical
  for (size_t j = 0; j < s.size() && j < maxVal.size(); j++)
  {
    size_t pidx = maxVal[j].size() > 1 ? LocalSystem::patch : 0;
    if (pidx < maxVal[j].size() && fabs(s[j]) > fabs(maxVal[j][pidx].second))
      maxVal[j][pidx] = std::make_pair(X,s[j]);
  }

  return true;
}


bool Elasticity::evalSol (Vector& s, const Vectors& eV, const FiniteElement& fe,
                          const Vec3& X, bool toLocal, Vec3* pdir) const
{
  if (eV.empty())
  {
    std::cerr <<" *** Elasticity::evalSol: No solutions vector."<< std::endl;
    return false;
  }
  else if (!eV.front().empty() && eV.front().size() != nsd*fe.dNdX.rows())
  {
    std::cerr <<" *** Elasticity::evalSol: Invalid displacement vector."
	      <<"\n     size(eV) = "<< eV.front().size() <<"   size(dNdX) = "
	      << fe.dNdX.rows() <<","<< fe.dNdX.cols() << std::endl;
    return false;
  }

  // Evaluate the deformation gradient, dUdX, and/or the strain tensor, eps
  Matrix Bmat;
  Tensor dUdX(nDF);
  SymmTensor eps(nsd,axiSymmetry);
  if (!this->kinematics(eV.front(),fe.N,fe.dNdX,X.x,Bmat,dUdX,eps))
    return false;

  // Add strains due to temperature expansion, if any
  double epsT = this->getThermalStrain(eV.back(),fe.N,X);
  if (epsT != 0.0) eps -= epsT;

  SymmTensor sigma(nsd, axiSymmetry || material->isPlaneStrain());
  if (wantStrain)
    sigma.copy(eps);
  else
  {
    // Calculate the stress tensor through the constitutive relation
    Matrix Cmat;
    double U = 0.0;
    if (!material->evaluate(Cmat,sigma,U,fe,X,dUdX,eps))
      return false;
    else if (epsT != 0.0 && nsd == 2 && material->isPlaneStrain())
      sigma(3,3) -= material->getStiffness(X)*epsT;
  }

  Vec3 p;
  bool havePval = false;
  if (toLocal && wantPrincipalStress)
  {
    // Calculate principal stresses and associated direction vectors
    if (sigma.size() == 4)
    {
      SymmTensor tmp(2); tmp = sigma; // discard the sigma_zz component
      havePval = pdir ? tmp.principal(p,pdir,2) : tmp.principal(p);
    }
    else
      havePval = pdir ? sigma.principal(p,pdir,2) : sigma.principal(p);

    // Congruence transformation to local coordinate system at current point
    if (locSys) sigma.transform(locSys->getTmat(X));
  }

  s = sigma;

  if (toLocal)
    s.push_back(sigma.vonMises());

  if (havePval)
  {
    s.push_back(p.x);
    s.push_back(p.y);
    if (sigma.dim() == 3)
      s.push_back(p.z);
  }

  return true;
}


bool Elasticity::evalSol (Vector& s, const STensorFunc& asol,
			  const Vec3& X) const
{
  s = asol(X);
  SymmTensor sigma(s);
  s = sigma;
  s.push_back(sigma.vonMises());

  if (wantPrincipalStress)
  {
    Vec3 p;
    sigma.principal(p);
    s.resize(s.size()+material->getNoIntVariables());
    s.push_back(p.x);
    s.push_back(p.y);
    if (nsd == 3)
      s.push_back(p.z);
  }
  return true;
}


bool Elasticity::evalEps (Vector& s, const Vector& eV, const FiniteElement& fe,
                          const Vec3& X) const
{
  if (eV.size() != fe.dNdX.rows()*nsd)
  {
    std::cerr <<" *** Elasticity::evalEps: Invalid displacement vector."
	      <<"\n     size(eV) = "<< eV.size() <<"   size(dNdX) = "
	      << fe.dNdX.rows() <<","<< fe.dNdX.cols() << std::endl;
    return false;
  }

  // Evaluate the strain tensor
  Matrix Bmat;
  SymmTensor eps(nsd,axiSymmetry);
  if (!this->kinematics(eV,fe.N,fe.dNdX,X.x,Bmat,eps,eps))
    return false;

  s = eps;
  return true;
}


Vector* Elasticity::getExtractionField (size_t ifield)
{
  return dS && ifield < primsol.size() ? &primsol[ifield] : nullptr;
}


bool Elasticity::getPrincipalDir (Matrix& pdir, size_t nPt, size_t idx) const
{
  if (!pDirBuf || idx < 1 || idx > 2) return false;

  if (pDirBuf->size() != nPt*2)
  {
    std::cerr <<" *** Elasticity::getPrincipalDir: Result point mismatch, nPt="
              << nPt <<", pDirBuf->size()="<< pDirBuf->size() << std::endl;
    return false;
  }

  pdir.resize(nsd,nPt);
  for (size_t i = 0; i < nPt; i++)
    pdir.fillColumn(1+i,(*pDirBuf)[2*i+idx-1].ptr());

  return true;
}


size_t Elasticity::getNoFields (int fld) const
{
  if (fld < 2)
    return nsd; // Displacement components

  size_t nf = nsd*(nsd+1)/2; // Symmetric stress tensor components
  if (nsd == 2 && (axiSymmetry || material->isPlaneStrain()))
    ++nf; // Include Hoop or normal stress

  if (fld == 2)
  {
    // Include von Mises stress and internal variables
    nf += 1 + material->getNoIntVariables();
    if (wantPrincipalStress)
      nf += nsd; // Include principal stress components
  }

#if INT_DEBUG > 1
  std::cout <<"Elasticity::getNoFields: "<< nf << std::endl;
#endif
  return nf;
}


std::string Elasticity::getField1Name (size_t i, const char* prefix) const
{
  if (i > nsd) i = 4;

  static const char* s[5] = { "u_x", "u_y", "u_z", "u_r", "displacement" };
  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[axiSymmetry ? 3-i : i];
}


std::string Elasticity::getField2Name (size_t i, const char* prefix) const
{
  size_t nVars = this->getNoFields(2);
  if (i >= nVars) return "";

  static const char* r[4] = { "s_rr", "s_zz", "s_tt", "s_zr" };
  static const char* s[6] = { "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz" };

  std::string name;
  if (prefix)
    name = std::string(prefix) + " ";
  else
    name.clear();

  // Number of components in the stress vector of this problem
  size_t nStress = this->getNoFields(3);

  if (nsd == 1)
    name += "Axial stress";
  else if (i == 2 && nStress == 3)
    name += s[3]; // No s_zz when plane stress
  else if (i < nStress)
  {
    size_t j = name.size();
    name += axiSymmetry ? r[i] : s[i];
    if (wantStrain) name[j] = 'e';
  }
  else if (i == nStress)
    name += "von Mises " + std::string(wantStrain ? "strain" : "stress");
  else if ((int)(i -= nStress) <= material->getNoIntVariables())
  {
    char varName[32];
    material->getInternalVariable(i,varName);
    name += varName;
  }
  else if ((i -= material->getNoIntVariables()) <= nsd)
    name += "P" + std::to_string(i);

  return name;
}


void Elasticity::initMaxVals (size_t nP)
{
#ifdef INT_DEBUG
  std::cout <<"Elasticity::initMaxVals: "<< maxVal.size()
            <<" --> "<< nP << std::endl;
#endif
  if (maxVal.empty() && nP > 0)
    maxVal.resize(this->getNoFields(2),PointValues(nP,PointValue(Vec3(),0.0)));
  else for (PointValues& pval : maxVal)
    std::fill(pval.begin(),pval.end(),PointValue(Vec3(),0.0));
}


void Elasticity::printMaxVals (std::streamsize precision, size_t comp) const
{
  size_t i1 = 1, i2 = maxVal.size();
  if (comp > i2)
    return;
  else if (comp > 0)
    i1 = i2 = comp;

  std::streamsize fldWidth = 8 + precision;
  utl::LogStream& os = IFEM::cout;
  for (size_t i = i1-1; i < i2; i++)
  {
    if (maxVal[i].empty())
      continue;
    else if (maxVal[i].size() == 1 && maxVal[i].front().second == 0.0)
      continue; // no value

    std::string name = this->getField2Name(i,nullptr);
    os <<"  Max "<< name <<":";
    if (name.size() < 16)
      os << std::string(16-name.size(),' ');
    for (size_t p = 0; p < maxVal[i].size(); p++)
    {
      std::streamsize oldPrec = os.precision(precision);
      std::ios::fmtflags oldF = os.flags(std::ios::scientific|std::ios::right);
      if (p > 0) os << std::string(23,' ');
      os << std::setw(fldWidth) << maxVal[i][p].second;
      os.precision(oldPrec);
      os.flags(oldF);
      os <<"  X = "<< maxVal[i][p].first;
      if (maxVal[i].size() > 1)
        os <<"\tPatch "<< p+1;
      os << std::endl;
    }
  }
}


NormBase* Elasticity::getNormIntegrand (AnaSol* asol) const
{
  if (asol && asol->getStressSol())
    return new ElasticityNorm(*const_cast<Elasticity*>(this),
                              asol->getStressSol(), asolProject ? 3 : 2);
  else
    return new ElasticityNorm(*const_cast<Elasticity*>(this));
}


ForceBase* Elasticity::getForceIntegrand (const Vec3* X0, AnaSol* asol) const
{
   return new ElasticityForce(*const_cast<Elasticity*>(this),X0,asol);
}


ForceBase* Elasticity::getForceIntegrand () const
{
   return new ElasticityForce(*const_cast<Elasticity*>(this));
}


size_t Elasticity::getVolumeIndex (bool withAnaSol) const
{
  return withAnaSol ? 5 : 3;
}


ElasticityNorm::ElasticityNorm (Elasticity& p, STensorFunc* a, int fld)
  : NormBase(p), anasol(a)
{
  nrcmp = myProblem.getNoFields(fld);
}


bool ElasticityNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
			      const Vec3& X) const
{
  Elasticity& problem = static_cast<Elasticity&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the inverse constitutive matrix at this point
  Matrix Cinv;
  if (!problem.formCmat(Cinv,fe,X,!Elasticity::wantStrain))
    return false;

  // Evaluate the finite element stress field
  Vector sigmah, sigma, error;
  if (!problem.evalSol(sigmah,pnorm.vec,fe,X))
    return false;

  // Evaluate the dual strain fields
  Vectors epsz(problem.numExtrFunction());
  for (size_t i = 1; i < pnorm.vec.size() && i-1 < epsz.size(); i++)
    if (!pnorm.vec[i].empty())
      if (!problem.evalEps(epsz[i-1],pnorm.vec[i],fe,X))
        return false;

  bool planeStrain = sigmah.size() == 4 && Cinv.rows() == 3;
  if (planeStrain)
  {
    // Remove the sigma_zz
    sigmah(3) = sigmah(4);
    sigmah.resize(3,utl::RETAIN);
  }

  double detJW = fe.detJxW;
  if (problem.isAxiSymmetric())
    detJW *= 2.0*M_PI*X.x;

  size_t ip = 0;
  // Integrate the energy norm a(u^h,u^h)
  pnorm[ip++] += sigmah.dot(Cinv*sigmah)*detJW;

  if (problem.haveLoads())
  {
    // Evaluate the body load
    Vec3 f = problem.getBodyforce(X);
    // Evaluate the displacement field
    Vec3 u = problem.evalSol(pnorm.vec.front(),fe.N);
    // Integrate the external energy (f,u^h)
    pnorm[ip] += f*u*detJW;
  }
  ip++;

  if (anasol)
  {
    // Evaluate the analytical stress field
    sigma = (*anasol)(X);
    if (sigma.size() == 4 && Cinv.rows() == 3)
    {
      // Remove the sigma_zz when plane strain
      sigma(3) = sigma(4);
      sigma.resize(3,utl::RETAIN);
    }

    // Integrate the energy norm a(u,u)
    pnorm[ip++] += sigma.dot(Cinv*sigma)*detJW;
    // Integrate the error in energy norm a(u-u^h,u-u^h)
    error = sigma - sigmah;
    pnorm[ip++] += error.dot(Cinv*error)*detJW;
  }

  // Integrate the volume
  pnorm[ip++] += detJW;

  for (const Vector& eps : epsz)
  {
    // Evaluate the variational-consistent postprocessing quantity, a(u^h,w)
    // (typically a sectional force component, or a mean stress component)
    pnorm[ip++] += sigmah.dot(eps)*detJW;
    if (anasol) // Evaluate the corresponding exact quantity, a(u,w)
      pnorm[ip++] += sigma.dot(eps)*detJW;
  }

#if INT_DEBUG > 3
  std::cout <<"\nElasticityNorm::evalInt(X = "<< X <<")\nsigma^h ="<< sigmah;
  if (anasol) std::cout <<"sigma ="<< sigma;
  size_t jp = anasol ? 5 : 3;
  for (size_t i = 0; i < epsz.size(); i++, jp += (anasol ? 2 : 1))
    std::cout <<"epsz("<< i+1 <<") ="<< epsz[i]
              <<"a(u^,w"<< i+1 <<"): "<< pnorm[jp] << std::endl;
#endif

  size_t j, k;
  for (const Vector& psol : pnorm.psol)
    if (!psol.empty())
    {
      // Evaluate projected stress field
      Vector sigmar(sigmah.size());
      for (j = k = 0; j < nrcmp && k < sigmar.size(); j++)
	if (!planeStrain || j != 2)
	  sigmar[k++] = psol.dot(fe.N,j,nrcmp);

      // Integrate the energy norm a(u^r,u^r)
      pnorm[ip++] += sigmar.dot(Cinv*sigmar)*detJW;
      // Integrate the error in energy norm a(u^r-u^h,u^r-u^h)
      error = sigmar - sigmah;
      pnorm[ip++] += error.dot(Cinv*error)*detJW;

      if (Elasticity::wantStrain)
      {
        // Convert to stresses, for L2-enorm
        sigmar = Cinv*sigmar;
        error  = Cinv*error;
      }
      double l2u = sigmar.norm2();
      double l2e = error.norm2();

      // Integrate the L2-norm (sigma^r,sigma^r)
      pnorm[ip++] += l2u*l2u*detJW;
      // Integrate the error in L2-norm (sigma^r-sigma^h,sigma^r-sigma^h)
      pnorm[ip++] += l2e*l2e*detJW;

      if (anasol)
      {
	// Integrate the error in the projected solution a(u-u^r,u-u^r)
	error = sigma - sigmar;
	pnorm[ip++] += error.dot(Cinv*error)*detJW;
	ip++; // Make room for the local effectivity index here
      }
    }

  return true;
}


bool ElasticityNorm::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
			      const Vec3& X, const Vec3& normal) const
{
  const Elasticity& problem = static_cast<const Elasticity&>(myProblem);
  if (!problem.haveLoads()) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the surface traction
  Vec3 T = problem.getTraction(X,normal);
  // Evaluate the displacement field
  Vec3 u = problem.evalSol(pnorm.vec.front(),fe.N);

  double detJW = fe.detJxW;
  if (problem.isAxiSymmetric())
    detJW *= 2.0*M_PI*X.x;

  // Integrate the external energy
  pnorm[1] += T*u*detJW;
  return true;
}


bool ElasticityNorm::finalizeElement (LocalIntegral& elmInt)
{
  if (!anasol) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate local effectivity indices as a(e^r,e^r)/a(e,e)
  // with e^r = u^r - u^h  and  e = u - u^h
  for (size_t ip = this->getNoFields(1); ip+5 < pnorm.size(); ip += 6)
    pnorm[ip+5] = pnorm[ip+1] / pnorm[3];

  return true;
}


size_t ElasticityNorm::getNoFields (int group) const
{
  if (group == 0)
    return this->NormBase::getNoFields();
  else if (group == 1 || group == -1)
  {
    size_t nExt = static_cast<Elasticity&>(myProblem).numExtrFunction();
    return anasol ? 5 + 2*nExt : 3 + nExt;
  }
  else if (group > 0 || !prjsol[-group-2].empty())
    return anasol ? 6 : 4;
  else
    return 0;
}


std::string ElasticityNorm::getName (size_t i, size_t j,
                                     const char* prefix) const
{
  size_t nx = i > 1 ? 1 : static_cast<Elasticity&>(myProblem).numExtrFunction();
  if (i == 1 && anasol) nx *= 2;
  if (i == 0 || j == 0 || j > 5+nx)
    return this->NormBase::getName(i,j,prefix);

  static const char* u[9] = {
    "a(u^h,u^h)^0.5",
    "((f,u^h)+(t,u^h))^0.5",
    "a(u,u)^0.5",
    "a(e,e)^0.5, e=u-u^h",
    "volume",
    "a(u^h,w)",
    "a(u,w)",
    "a(z^h,z^h)^0.5",
    "(E(u)*E(z))^0.5, E(v)=a(e,e), e=v^r-v^h"
  };

  static const char* p[6] = {
    "a(u^r,u^r)^0.5",
    "a(e,e)^0.5, e=u^r-u^h",
    "(u^r,u^r)^0.5",
    "(e,e)^0.5, e=u^r-u^h",
    "a(e,e)^0.5, e=u-u^r",
    "effectivity index"
  };

  const char** s = i > 1 ? p : u;
  if (i == 1)
  {
    if (anasol)
      nx /= 2;
    else if (j > 2 && j < 4)
      j += 2;
    if (j > 5 && nx > 1)
    {
      j -= 5;
      char comp[32];
      if (!anasol)
        sprintf(comp,"a(u^h,w%zu)",j);
      else if (j%2)
        sprintf(comp,"a(u^h,w%zu)",j/2+1);
      else
        sprintf(comp,"a(u,w%zu)",j/2);
      return comp;
    }
    else if (j <= 2 && prefix)
      if (!strncmp(prefix,"Dual",4))
        return s[6+j];
  }
  if (!prefix)
    return s[j-1];

  return prefix + std::string(" ") + s[j-1];
}


bool ElasticityNorm::hasElementContributions (size_t i, size_t j) const
{
  return i > 1 || j != 2;
}


LocalIntegral* ElasticityForce::getLocalIntegral (size_t nen, size_t iEl,
                                                  bool) const
{
  if (nodal) {
    ElmMats* result = new ElmMats(false);
    result->resize(0,1);
    result->redim(static_cast<Elasticity&>(myProblem).getNoSpaceDim()*nen);
    return result;
  }

  return this->ForceBase::getLocalIntegral(nen,iEl);
}


bool ElasticityForce::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
			       const Vec3& X, const Vec3& normal) const
{
  const Elasticity& problem = static_cast<const Elasticity&>(myProblem);
  size_t nsd = fe.dNdX.cols();

  // Numerical approximation of stress
  Vector stress;
  problem.evalSol(stress,elmInt.vec,fe,X);
  SymmTensor sigmah(stress);

  // Finite element traction
  Vec3 th = sigmah*normal;

  if (nodal) // Evaluate nodal force vectors
    return this->evalForce(static_cast<ElmMats&>(elmInt),th,fe);

  // Evaluate integrated forces (and errors) over the element
  return this->evalForce(static_cast<ElmNorm&>(elmInt),
                         th,X,normal,fe.detJxW,nsd);
}


bool ElasticityForce::evalForce (ElmNorm& pnorm, const Vec3& th,
				 const Vec3& X, const Vec3& normal,
				 double detJW, size_t nsd) const
{
  // Numerical force term
  size_t i, ip = 0;
  for (i = 0; i < nsd; i++)
    pnorm[ip++] += th[i]*detJW;

  if (X0) {
    // Numerical torque term
    Vec3 T(X-(*X0),th);
    if (nsd == 2)
      pnorm[ip] += T[2]*detJW;
    else for (i = 0; i < nsd; i++)
      pnorm[ip++] += T[i]*detJW;
  }

  return true;
}


bool ElasticityForce::evalForce (ElmMats& elmat, const Vec3& th,
                                 const FiniteElement& fe) const
{
  // Integrate the nodal force vector
  Vector& ES = elmat.b.front();
  size_t nsd = fe.dNdX.cols();
  for (size_t a = 1; a <= fe.N.size(); a++)
    for (size_t d = 1; d <= nsd; d++)
      ES(nsd*(a-1)+d) += th[d-1]*fe.N(a)*fe.detJxW;

  return true;
}


size_t ElasticityForce::getNoComps () const
{
  size_t nsd = static_cast<Elasticity&>(myProblem).getNoSpaceDim();
  if (nodal) return nsd;

  size_t nf = nsd; // traction components

  if (X0)
    nf += nsd == 2 ? 1 : nsd; // add torque component(s)

  return nf;
}
