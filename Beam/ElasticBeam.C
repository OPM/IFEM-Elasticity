// $Id$
//==============================================================================
//!
//! \file ElasticBeam.C
//!
//! \date Aug 10 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class representing a nonlinear elastic beam.
//!
//==============================================================================

#include "ElasticBeam.h"
#include "FiniteElement.h"
#include "HHTMats.h"
#include "ElmNorm.h"
#include "AnaSol.h"
#include "Functions.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include "tinyxml.h"


/*!
  \brief Base class for beam FE problems containing data for local element axes.
*/

class LocalBeamAxes
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  LocalBeamAxes() : Tlg(3,3) {}

public:
  //! \brief Empty destructor.
  virtual ~LocalBeamAxes() {}

  Matrix Tlg; //!< Local to global transformation matrix
};


/*!
  \brief Class collecting element matrices and data for a beam FE problem.

  \details This class is defined as a template, where the template argument
  defines the parent class. This way sub-classes of different parents may be
  generated compile-time from the same source code.
*/

template <class T> class BeamMats : public T, public LocalBeamAxes
{
public:
  //! \brief Default constructor.
  BeamMats(double a1 = 0.0, double a2 = 0.0, double b = 0.0, double c = 0.0) {}
  //! \brief Empty destructor.
  virtual ~BeamMats() {}
};

//! Element matrices and data for static beam FE problems
typedef BeamMats<ElmMats> BeamElmMats;

//! Element matrices and data for linear dynamic beam FE problems
typedef BeamMats<NewmarkMats> NewmarkBeamMats;

template<> NewmarkBeamMats::BeamMats (double a1, double a2, double b, double c)
  : NewmarkMats(a1,a2,b,c) {}

//! Element matrices and data for nonlinear dynamic beam FE problems
typedef BeamMats<HHTMats> HHTBeamMats;

template<> HHTBeamMats::BeamMats (double a1, double a2, double b, double)
  : HHTMats(a1,a2,b,true) {}


ElasticBeam::ElasticBeam (unsigned short int n) : inLocalAxes(true)
{
  npv = 6; // Number of primary unknowns per node
  nSV = n; // Number of solution vectors in core

  // Default material parameters
  E   = 2.05e11;
  G   = 7.94e10;
  rho = 7.85e3;

  this->initPropFunc();
}


ElasticBeam::~ElasticBeam ()
{
  delete lineLoad;
  delete cplLoad;
}


void ElasticBeam::initPropFunc ()
{
  myProp = nullptr;
  lineLoad = cplLoad = nullptr;
}


void ElasticBeam::printLog () const
{
  IFEM::cout <<"ElasticBeam: E = "<< E <<", G = "<< G <<", rho = "<< rho;
  if (!myProp)
  {
    // Using default cross section properties
    static BeamProperty defProp;
    const_cast<ElasticBeam*>(this)->myProp = &defProp;
  }

  IFEM::cout <<"\n             A = "<< myProp->A
             <<", Ix = "<< myProp->Ix
             <<", Iy = "<< myProp->Iy
             <<", Iz = "<< myProp->Iz <<", It = "<< myProp->It
             <<"\n             Ky = "<< myProp->Ky
             <<", Kz = "<< myProp->Kz
             <<", Sy = "<< myProp->Sy <<", Sz = "<< myProp->Sz
             << std::endl;
}


void ElasticBeam::setBeamLoad (VecFunc* load, bool cpl)
{
  if (cpl)
    cplLoad = load;
  else
    lineLoad = load;
}


void ElasticBeam::parseBeamLoad (const TiXmlElement* load)
{
  std::string type;
  utl::getAttribute(load,"type",type,true);
  if (!strcasecmp(load->Value(),"lineload"))
  {
    IFEM::cout <<"\tDistributed load";
    lineLoad = utl::parseVecFunc(load->FirstChild()->Value(),type);
  }
  else if (!strcasecmp(load->Value(),"cplload"))
  {
    IFEM::cout <<"\tDistributed moment load";
    cplLoad = utl::parseVecFunc(load->FirstChild()->Value(),type);
  }
  else
    return;

  IFEM::cout << std::endl;
}


LocalIntegral* ElasticBeam::getLocalIntegral (size_t, size_t iEl, bool) const
{
  ElmMats* result = nullptr;
  if (this->inActive(iEl))
    return result; // element is not in current material group

  if (m_mode != SIM::DYNAMIC)
    result = new BeamElmMats();
  else if (intPrm[3] > 0.0)
    result = new NewmarkBeamMats(intPrm[0],intPrm[1],intPrm[2],intPrm[3]);
  else
    result = new HHTBeamMats(intPrm[2],intPrm[0],intPrm[1]);

  switch (m_mode)
  {
    case SIM::STATIC:
    case SIM::MASS_ONLY:
      result->resize(1,1);
      break;

    case SIM::ARCLEN:
      result->resize(1,2);
      break;

    case SIM::DYNAMIC:
      result->resize(intPrm[3] >= 0.0 ? 3 : 4, intPrm[3] > 0.0 ? 1 : 2);
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

  result->redim(12);
  return result;
}


bool ElasticBeam::initElement (const std::vector<int>& MNPC,
                               const FiniteElement& fe, const Vec3&, size_t,
                               LocalIntegral& elmInt)
{
  if (!this->initElement(MNPC,elmInt))
    return false;

  Vec3 e1 = fe.XC[1] - fe.XC[0]; // Initial local X-axis
  const Vector& eV = elmInt.vec.front();
  if (!eV.empty())
  {
    // Fetch nodal displacements
    Vec3 U0(eV.ptr());
    Vec3 U1(eV.ptr()+eV.size()-npv);
    e1 += U1 - U0; // Deformed local X-axis
  }

  // Calculate the co-rotated element coordinate system
  if (e1.normalize() <= 1.0e-8)
  {
    std::cerr <<" *** ElasticBeam::initElement: Zero beam length"<< std::endl;
    return false;
  }

  if (fe.Tn.size() < 2)
  {
    std::cerr <<" *** ElasticBeam::initElement: No end rotations"<< std::endl;
    return false;
  }

  Vec3 e2(fe.Tn[0][1]+fe.Tn[1][1]); // Sum of the nodal Y-axes
  Vec3 e3(fe.Tn[0][2]+fe.Tn[1][2]); // Sum of the nodal Z-axes
  if (fabs(e3.x) < fabs(e2.x))
  {
    // Local Y-axis = e3xe1 / |e3xe1|
    e2.cross(fe.Te*e3,e1);
    e2.normalize();
    // Local Z-axis = e1xe2
    e3.cross(e1,e2);
  }
  else
  {
    // Local Z-axis = e1xe2 / |e1xe2|
    e3.cross(e1,fe.Te*e2);
    e3.normalize();
    // Local Y-axis = e3xe1
    e2.cross(e3,e1);
  }

  Matrix& Tlg = this->getLocalAxes(elmInt);
  Tlg.fillColumn(1,e1.ptr());
  Tlg.fillColumn(2,e2.ptr());
  Tlg.fillColumn(3,e3.ptr());

#if INT_DEBUG > 1
  std::cout <<"ElasticBeam: local-to-global transformation matrix:"<< Tlg
            <<"ElasticBeam: T1n\n"<< fe.Tn[0] <<"ElasticBeam: T2n\n"<< fe.Tn[1];
#endif

  return true;
}


/*!
  This method is just a workaround for that we can not cast directly from
  LocalIntegral& to LocalBeamAxes& due to that class LocalBeamAxes does not
  inherit class LocalIntegral.
*/

Matrix& ElasticBeam::getLocalAxes (LocalIntegral& elmInt) const
{
  if (m_mode != SIM::DYNAMIC)
    return static_cast<BeamElmMats&>(elmInt).Tlg;
  else if (intPrm[3] > 0.0)
    return static_cast<NewmarkBeamMats&>(elmInt).Tlg;
  else
    return static_cast<HHTBeamMats&>(elmInt).Tlg;
}


/*!
  \details This matrix is taken from page 12 in
  http://people.duke.edu/~hpgavin/cee421/frame-finite-def.pdf
*/

void ElasticBeam::getMaterialStiffness (Matrix& EK, double L,
                                        double EA,  double GIt,
                                        double EIy, double EIz,
                                        double Aly, double Alz) const
{
  double L3 = L*L*L;

  EK.resize(12,12,true);
  EK( 1, 1) =  EK( 7, 7) =  EA/L;
  EK( 2, 2) =  EK( 8, 8) =  12.0*EIy/(L3+L3*Aly);
  EK( 3, 3) =  EK( 9, 9) =  12.0*EIz/(L3+L3*Alz);
  EK( 4, 4) =  EK(10,10) =  GIt/L;
  EK( 5, 5) =  EK(11,11) =  EIz*(4.0+Alz)/(L+L*Alz);
  EK( 6, 6) =  EK(12,12) =  EIy*(4.0+Aly)/(L+L*Aly);

  EK( 3, 5) =  EK( 3,11) = -EK(3,3)*L*0.5;
  EK( 2, 6) =  EK( 2,12) =  EK(2,2)*L*0.5;
  EK( 1, 7) = -EK( 1, 1);

  EK( 6, 8) =  EK( 8,12) = -EK(2,6);
  EK( 5, 9) =  EK( 9,11) = -EK(3,5);
  EK( 4,10) = -EK( 4, 4);

  EK( 2, 8) = -EK( 2, 2);
  EK( 3, 9) = -EK( 3, 3);
  EK( 5,11) =  EIz*(2.0-Alz)/(L+L*Alz);
  EK( 6,12) =  EIy*(2.0-Aly)/(L+L*Aly);

  // Lower triangle from symmetry
  size_t i, j;
  for (i = 2; i <= 12; i++)
    for (j = 1; j < i; j++)
      EK(i,j) = EK(j,i);

  if (fabs(myProp->Sy) + fabs(myProp->Sz) > 0.00001*sqrt(myProp->A))
  {
    // Adjustment due to non-symmetric cross section
    for (j = 1; j <= 12; j++)
      for (i = 0; i <= 6; i += 6)
        EK(i+4,j) -= myProp->Sz*EK(i+2,j) + myProp->Sy*EK(i+3,j);

    for (i = 1; i <= 12; i++)
      for (j = 0; j <= 6; j += 6)
        EK(i,j+4) -= EK(i,j+2)*myProp->Sz + EK(i,j+3)*myProp->Sy;
  }

#if INT_DEBUG > 1
  std::cout <<"\nElasticBeam: local material stiffness matrix:"<< EK;
#endif
}


/*!
  \details This matrix is taken from page 14 in
  http://people.duke.edu/~hpgavin/cee421/frame-finite-def.pdf
*/

void ElasticBeam::getGeometricStiffness (Matrix& EK, double N, double L,
                                         double EIy, double EIz,
                                         double Aly, double Alz) const
{
  double L2 = L*L;
  double Cy = N/(L*(1.0+Aly)*(1.0+Aly));
  double Cz = N/(L*(1.0+Alz)*(1.0+Alz));

  EK.resize(12,12,true);
  EK( 2, 2) =  EK( 8, 8) =  Cy*(1.2+Aly*(2.0+Aly));
  EK( 3, 3) =  EK( 9, 9) =  Cz*(1.2+Alz*(2.0+Alz));
  EK( 4, 4) =  EK(10,10) =  N*myProp->It/(myProp->A*L);
  EK( 5, 5) =  EK(11,11) =  Cz*L2*(0.4+Alz*(0.5+0.25*Alz))/3.0;
  EK( 6, 6) =  EK(12,12) =  Cy*L2*(0.4+Aly*(0.5+0.25*Aly))/3.0;

  EK( 3, 5) =  EK( 3,11) = -Cz*L*0.1;
  EK( 5, 9) =  EK( 9,11) = -EK(3,5);
  EK( 2, 6) =  EK( 2,12) =  Cy*L*0.1;
  EK( 6, 8) =  EK( 8,12) = -EK(2,6);

  EK( 2, 8) = -EK( 2, 2);
  EK( 3, 9) = -EK( 3, 3);
  EK( 4,10) = -EK( 4, 4);

  EK( 5,11) = -Cz*L2*(0.1+Alz*(0.5+0.25*Alz))/3.0;
  EK( 6,12) = -Cy*L2*(0.1+Aly*(0.5+0.25*Aly))/3.0;

  // Lower triangle from symmetry
  size_t i, j;
  for (i = 2; i <= 12; i++)
    for (j = 1; j < i; j++)
      EK(i,j) = EK(j,i);

  if (fabs(myProp->Sy) + fabs(myProp->Sz) > 0.00001*sqrt(myProp->A))
  {
    // Adjustment due to non-symmetric cross section
    for (j = 1; j <= 12; j++)
      for (i = 0; i <= 6; i += 6)
        EK(i+4,j) -= myProp->Sz*EK(i+2,j) + myProp->Sy*EK(i+3,j);

    for (i = 1; i <= 12; i++)
      for (j = 0; j <= 6; j += 6)
        EK(i,j+4) -= EK(i,j+2)*myProp->Sz + EK(i,j+3)*myProp->Sy;
  }

#if INT_DEBUG > 1
  std::cout <<"ElasticBeam: local geometric stiffness matrix:"<< EK;
#endif
}


/*!
  \details This matrix is taken from
  http://www.stanford.edu/class/aa244a/mfiles/mbeam3d.m
*/

void ElasticBeam::getMassMatrix (Matrix& EM, double rhoA, double Ixx,
                                 double Iyy, double Izz, double L) const
{
  double AM  = L*rhoA;
  double AM2 = L*AM;
  double AM3 = L*AM2;

  EM.resize(12,12,true);
  EM( 1, 1) =  EM( 7, 7) = AM/3.0;
  EM( 2, 2) =  EM( 8, 8) = AM*1.3/3.5 + Izz*1.2/L;
  EM( 3, 3) =  EM( 9, 9) = AM*1.3/3.5 + Iyy*1.2/L;
  EM( 4, 4) =  EM(10,10) = L*Ixx/3.0;
  EM( 5, 5) =  EM(11,11) = AM3/105.0  + Iyy*L*0.4/3.0;
  EM( 6, 6) =  EM(12,12) = AM3/105.0  + Izz*L*0.4/3.0;

  EM( 1, 7) =  AM/6.0;
  EM( 2, 6) =  AM2*1.1/21.0 + Izz*0.1;
  EM( 2, 8) =  AM *0.9/7.0  - Izz*1.2/L;
  EM( 2,12) = -AM2*1.3/42.0 + Izz*0.1;
  EM( 3, 5) = -AM2*1.1/21.0 - Iyy*0.1;
  EM( 3, 9) =  AM *0.9/7.0  - Iyy*1.2/L;
  EM( 3,11) =  AM2*1.3/42.0 - Iyy*0.1;
  EM( 4,10) =  EM(4,4)*0.5;
  EM( 5, 9) = -EM(3,11);
  EM( 5,11) = -AM3/140.0    - Iyy*L/30.0;
  EM( 6, 8) = -EM(2,12);
  EM( 6,12) = -AM3/140.0    - Izz*L/30.0;
  EM( 8,12) = -EM(2,6);
  EM( 9,11) = -EM(3,5);

  // Lower triangle from symmetry
  size_t i, j;
  for (i = 2; i <= 12; i++)
    for (j = 1; j < i; j++)
      EM(i,j) = EM(j,i);

#if INT_DEBUG > 1
  std::cout <<"ElasticBeam: local mass matrix:"<< EM;
#endif
}


/*!
  \brief Eccentricity transformation of a beam element matrix.
*/

static void eccTransform (Matrix& A, const Vec3& e1, const Vec3& e2)
{
  for (int j = 1; j <= 12; j++)
  {
    A( 4,j) +=  e1.z*A(2,j) - e1.y*A(3,j);
    A( 5,j) += -e1.z*A(1,j) + e1.x*A(3,j);
    A( 6,j) +=  e1.y*A(1,j) - e1.x*A(2,j);
    A(10,j) +=  e2.z*A(8,j) - e2.y*A(9,j);
    A(11,j) += -e2.z*A(7,j) + e2.x*A(9,j);
    A(12,j) +=  e2.y*A(7,j) - e2.x*A(8,j);
  }

  for (int i = 1; i <= 12; i++)
  {
    A(i, 4) +=  e1.z*A(i,2) - e1.y*A(i,3);
    A(i, 5) += -e1.z*A(i,1) + e1.x*A(i,3);
    A(i, 6) +=  e1.y*A(i,1) - e1.x*A(i,2);
    A(i,10) +=  e2.z*A(i,8) - e2.y*A(i,9);
    A(i,11) += -e2.z*A(i,7) + e2.x*A(i,9);
    A(i,12) +=  e2.y*A(i,7) - e2.x*A(i,8);
  }
}


bool ElasticBeam::evalInt (LocalIntegral& elmInt,
                           const FiniteElement& fe,
                           const Vec3& X) const
{
  // Calculate initial element length
  Vec3 X0 = fe.XC[1] - fe.XC[0];
  double L0 = X0.length();
  if (L0 <= 1.0e-8)
  {
    std::cerr <<" *** ElasticBeam::evalInt: Zero initial element length "
              << L0 << std::endl;
    return false;
  }
#if INT_DEBUG > 1
  std::cout <<"ElasticBeam: iel = "<< fe.iel <<" X = "<< X <<", L0 = "<< L0;
#endif

  const Vector& eV = elmInt.vec.front();
  if (!eV.empty())
  {
    // Calculate current element length
    Vec3 U1(eV.ptr(),npv);
    Vec3 U2(eV.ptr()+eV.size()-npv,npv);
    Vec3 X1 = X0 + U2 - U1;
    double L = X1.length();
    if (L <= 1.0e-8)
    {
      std::cerr <<" *** ElasticBeam::evalInt: Zero element length "
                << L << std::endl;
      return false;
    }
#if INT_DEBUG > 1
    std::cout <<" L = "<< L <<", U1 = "<< U1 <<" U2 = "<< U2;
#endif
  }

  if (!myProp)
  {
    std::cerr <<" *** ElasticBeam::evalInt: No properties."<< std::endl;
    return false;
  }

  // Evaluate beam stiffness and mass properties at this point
  double EA, EIy, EIz, GIt, Aly, Alz;
  double rhoA, I_xx, I_yy, I_zz, CG_y, CG_z;
  bool hasGrF = gravity.isZero() ? false : eS > 0;
  myProp->eval(X, L0, E, G, rho, hasGrF, eM > 0,
               EA, EIy, EIz, GIt, Aly, Alz,
               rhoA, CG_y, CG_z, I_xx, I_yy, I_zz);
#if INT_DEBUG > 1
  std::cout <<"\n             EA = "<< EA
            <<" EI = "<< EIy <<" "<< EIz <<" GIt = "<< GIt
            <<", Alpha_y = "<< Aly <<" Alpha_z = "<< Alz;
  std::cout <<", rho*A = "<< rhoA <<" rho*I = "<< I_xx <<" "<< I_yy <<" "<< I_zz
            <<", CoG = "<< CG_y <<" "<< CG_z << std::endl;
#endif

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (hasGrF)
  {
    // External (gravitation) forces
    Vector& S = elMat.b[eS-1];
    Vec3 gvec = (0.5*rhoA*L0) * (gravity*this->getLocalAxes(elmInt));
    for (unsigned short int i = 1; i <= 3; i++)
      S(i) = S(npv+i) = gvec[i-1];

    // If the centre of gravity has an offset w.r.t. the neutral axis,
    // it will result in an additional torque load on the element
    if (CG_y != 0.0 || CG_z != 0.0)
      S(4) = S(npv+4) = gvec.z*CG_y - gvec.y*CG_z;
#if INT_DEBUG > 1
    std::cout <<"ElasticBeam: S_ext"<< S << std::endl;
#endif
  }

  if (eKm) // Evaluate the material stiffness matrix
    this->getMaterialStiffness(elMat.A[eKm-1],L0,EA,GIt,EIy,EIz,Aly,Alz);

  Vector v;
  double N = 0.0;

  if ((iS || eKg) && eV.normInf() > 1.0e-16*L0)
  {
    v = eV; // Transform the element displacement vector to local coordinates
    const Matrix& Tlg = this->getLocalAxes(elmInt);
    for (size_t k = 1; k < v.size(); k += 3)
      if (!utl::transform(v,Tlg,k,true))
        return false;

#if INT_DEBUG > 1
    std::cout <<"ElasticBeam: v"<< v;
#endif
    if (eKg) N = EA*(v(7)-v(1))/L0; // Axial force
  }

  if (iS && !v.empty())
  {
    v *= -1.0;

    // Internal forces, S_int = Km*v
    Matrix tmpKm;
    if (!eKm) this->getMaterialStiffness(tmpKm,L0,EA,GIt,EIy,EIz,Aly,Alz);
    Matrix& Km = eKm ? elMat.A[eKm-1] : tmpKm;
    if (!Km.multiply(v,elMat.b[iS-1],false,iS == eS))
      return false;

#if INT_DEBUG > 1
    if (iS == eS)
      std::cout <<"ElasticBeam: S_ext - S_int"<< elMat.b[iS-1] << std::endl;
    else
      std::cout <<"ElasticBeam: -S_int"<< elMat.b[iS-1] << std::endl;
#endif
  }

  if (eKg && N != 0.0)
  {
#if INT_DEBUG > 1
    std::cout <<"ElasticBeam: Axial force, N = "<< N << std::endl;
#endif

    // Evaluate the geometric stiffness matrix
    if (eKg == eKm)
    {
      Matrix Kg(12,12);
      this->getGeometricStiffness(Kg,N,L0,EIy,EIz,Aly,Alz);
      elMat.A[eKm-1].add(Kg);
    }
    else
      this->getGeometricStiffness(elMat.A[eKg-1],N,L0,EIy,EIz,Aly,Alz);
  }

  if (eM)
  {
    // Evaluate the mass matrix
    this->getMassMatrix(elMat.A[eM-1],rhoA,I_xx,I_yy,I_zz,L0);
    if (CG_y != 0.0 || CG_z != 0.0) // Transform to neutral axis location
      eccTransform(elMat.A[eM-1],Vec3(0.0,-CG_y,-CG_z),Vec3(0.0,-CG_y,-CG_z));
  }

  return true;
}


bool ElasticBeam::finalizeElement (LocalIntegral& elmInt,
                                   const TimeDomain& time, size_t iGP)
{
  if (inLocalAxes)
  {
    size_t i, k;
    ElmMats& elMat = static_cast<ElmMats&>(elmInt);
    const Matrix& Tlg = this->getLocalAxes(elmInt);

    // Transform the element matrices to global coordinates
    for (i = 0; i < elMat.A.size(); i++)
      for (k = 1; k < elMat.A[i].cols(); k += 3)
        if (!utl::transform(elMat.A[i],Tlg,k))
          return false;

    // Transform the element force vectors to global coordinates
    for (i = 0; i < elMat.b.size(); i++)
      for (k = 1; k < elMat.b[i].size(); k += 3)
        if (!utl::transform(elMat.b[i],Tlg,k))
          return false;
  }

  return this->ElasticBase::finalizeElement(elmInt,time,iGP);
}


Vec3 ElasticBeam::displacement (const FiniteElement& fe, const Vector& eV) const
{
  Vec3 d;
  for (size_t k = 0; k < 3; k++)
    d[k] = eV.dot(fe.N,k,npv);

  return d;
}


NormBase* ElasticBeam::getNormIntegrand (AnaSol* asol) const
{
  return new ElasticBeamNorm(*const_cast<ElasticBeam*>(this),
                             asol ? asol->getVectorSol() : nullptr);
}


size_t ElasticBeamNorm::getNoFields (int fld) const
{
  if (fld == 0)
    return 1;

  return anasol ? 3 : 1;
}


std::string ElasticBeamNorm::getName (size_t, size_t j, const char* prfix) const
{
  static const char* u[3] = {
     "||u^h||_L2",
     "||u||_L2",
     "||e||_L2, e=u^h-u"
  };

  if (!prfix)
    return u[j-1];

  return prfix + std::string(" ") + u[j-1];
}


bool ElasticBeamNorm::initElement (const std::vector<int>& MNPC,
                                   const FiniteElement&, const Vec3&, size_t,
                                   LocalIntegral& elmInt)
{
  return this->initProjection(MNPC,elmInt) &&
         myProblem.initElement(MNPC,elmInt);
}


bool ElasticBeamNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                               const Vec3& X) const
{
  const ElasticBeam& problem = static_cast<const ElasticBeam&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the FE displacements
  Vec3 Uh = problem.displacement(fe,elmInt.vec.front());

  // L2-norm of FE displacements
  pnorm[0] += Uh*Uh*fe.detJxW;

  if (anasol)
  {
    // Exact displacements
    Vec3 U = (*anasol)(X);

    // L2-norm of exact solution
    pnorm[1] += U*U*fe.detJxW;

    // L2-norm of displacement error
    U -= Uh;
    pnorm[2] += U*U*fe.detJxW;
  }

  return true;
}
