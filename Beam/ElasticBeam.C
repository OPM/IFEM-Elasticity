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
#include "BeamProperty.h"
#include "FiniteElement.h"
#include "HHTMats.h"
#include "ElmNorm.h"
#include "AnaSol.h"
#include "Functions.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include "tinyxml2.h"


namespace
{
  /*!
    \brief Base class for beam problems containing data for local element axes.
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

  template<class T> class BeamMats : public T, public LocalBeamAxes
  {
  public:
    //! \brief Default constructor.
    explicit BeamMats(double = 0.0, double = 0.0, double = 0.0, double = 0.0) {}
  };

  //! Element matrices and data for static beam FE problems
  using BeamElmMats = BeamMats<ElmMats>;
  //! Element matrices and data for linear dynamic beam FE problems
  using DynBeamMats = BeamMats<NewmarkMats>;
  //! Element matrices and data for nonlinear dynamic beam FE problems
  using HHTBeamMats = BeamMats<HHTMats>;

  //! \brief Instantiation for NewmarMats as parent class.
  template<> DynBeamMats::BeamMats (double a1, double a2, double b, double c)
    : NewmarkMats(a1,a2,b,c) {}

  //! \brief Instantiation for HHTMats as parent class.
  template<> HHTBeamMats::BeamMats (double a1, double a2, double b, double)
    : HHTMats(a1,a2,b,true) {}

} // namespace


ElasticBeam::ElasticBeam (unsigned short int n) : inLocalAxes(true)
{
  nsd = 3; // Number of spatial dimenstions
  npv = 6; // Number of primary unknowns per node
  nCS = n; // Number of consecutive solution states in core
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

  IFEM::cout << *myProp << std::endl;
}


void ElasticBeam::setBeamLoad (VecFunc* load, bool cpl)
{
  if (cpl)
    cplLoad = load;
  else
    lineLoad = load;
}


void ElasticBeam::parseBeamLoad (const tinyxml2::XMLElement* load)
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


BeamProperty* ElasticBeam::parseProp (const tinyxml2::XMLElement* prop)
{
  return new BeamProperty(prop);
}


int ElasticBeam::getIntegrandType() const
{
  int linearItgType = NO_DERIVATIVES | ELEMENT_CORNERS;
  return nSV > 0 ? linearItgType | NODAL_ROTATIONS : linearItgType;
}


LocalIntegral* ElasticBeam::getLocalIntegral (size_t, size_t iEl, bool) const
{
  ElmMats* result = nullptr;
  if (this->inActive(iEl))
    return result; // element is not in current material group

  if (m_mode != SIM::DYNAMIC)
    result = new BeamElmMats();
  else if (intPrm[3] > 0.0)
    result = new DynBeamMats(intPrm[0],intPrm[1],intPrm[2],intPrm[3]);
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
  result->Aname = matNames;
  result->Bname = vecNames;
  return result;
}


bool ElasticBeam::initElement (const std::vector<int>& MNPC,
                               const FiniteElement& fe, const Vec3&, size_t,
                               LocalIntegral& elmInt)
{
  if (!this->initElement(MNPC,elmInt))
    return false;

  Matrix& Tlg = this->getLocalAxes(elmInt);
  if (fe.Tn.size() < 2)
  {
    // No end rotation tensors - assuming linear analysis
    Tlg = fe.Te;
#if INT_DEBUG > 1
    std::cout <<"ElasticBeam: local-to-global transformation matrix:"<< Tlg;
#endif
    return true;
  }

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
  This method is just a workaround for that we can not cast directly from a
  LocalIntegral reference to a LocalBeamAxes reference, due to that class
  LocalBeamAxes does not inherit class LocalIntegral.
*/

Matrix& ElasticBeam::getLocalAxes (LocalIntegral& elmInt) const
{
  if (m_mode != SIM::DYNAMIC)
    return static_cast<BeamElmMats&>(elmInt).Tlg;
  else if (intPrm[3] > 0.0)
    return static_cast<DynBeamMats&>(elmInt).Tlg;
  else
    return static_cast<HHTBeamMats&>(elmInt).Tlg;
}


namespace // Private matrix transformation methods.
{
  //! \brief Eccentricity transformation of a beam element vector.
  void eccTransform (Vector& b, const Vec3& e1, const Vec3& e2)
  {
    b( 4) +=  e1.z*b(2) - e1.y*b(3);
    b( 5) += -e1.z*b(1) + e1.x*b(3);
    b( 6) +=  e1.y*b(1) - e1.x*b(2);
    b(10) +=  e2.z*b(8) - e2.y*b(9);
    b(11) += -e2.z*b(7) + e2.x*b(9);
    b(12) +=  e2.y*b(7) - e2.x*b(8);
  }

  //! \brief Eccentricity transformation of a beam element matrix.
  void eccTransform (Matrix& A, const Vec3& e1, const Vec3& e2)
  {
    for (size_t j = 1; j <= 12; j++)
    {
      A( 4,j) +=  e1.z*A(2,j) - e1.y*A(3,j);
      A( 5,j) += -e1.z*A(1,j) + e1.x*A(3,j);
      A( 6,j) +=  e1.y*A(1,j) - e1.x*A(2,j);
      A(10,j) +=  e2.z*A(8,j) - e2.y*A(9,j);
      A(11,j) += -e2.z*A(7,j) + e2.x*A(9,j);
      A(12,j) +=  e2.y*A(7,j) - e2.x*A(8,j);
    }

    for (size_t i = 1; i <= 12; i++)
    {
      A(i, 4) +=  e1.z*A(i,2) - e1.y*A(i,3);
      A(i, 5) += -e1.z*A(i,1) + e1.x*A(i,3);
      A(i, 6) +=  e1.y*A(i,1) - e1.x*A(i,2);
      A(i,10) +=  e2.z*A(i,8) - e2.y*A(i,9);
      A(i,11) += -e2.z*A(i,7) + e2.x*A(i,9);
      A(i,12) +=  e2.y*A(i,7) - e2.x*A(i,8);
    }
  }

  //! \brief Shear-centre transformation of a beam element matrix.
  void shearTransform (Matrix& A, double Sy, double Sz)
  {
    for (size_t j = 1; j <= 12; j++)
      for (size_t i = 0; i <= 6; i += 6)
        A(i+4,j) += Sy*A(i+3,j) - Sz*A(i+2,j);

    for (size_t i = 1; i <= 12; i++)
      for (size_t j = 0; j <= 6; j += 6)
        A(i,j+4) += A(i,j+3)*Sy - A(i,j+2)*Sz;
  }
}


/*!
  \details This matrix is taken from page 12 in
  http://people.duke.edu/~hpgavin/cee421/frame-finite-def.pdf
*/

void ElasticBeam::getMaterialStiffness (Matrix& EK, double L,
                                        double EA,  double GIt,
                                        double EIy, double EIz,
                                        double Aly, double Alz,
                                        double Sy,  double Sz) const
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
  for (size_t i = 2; i <= 12; i++)
    for (size_t j = 1; j < i; j++)
      EK(i,j) = EK(j,i);

  // Adjustment due to non-symmetric cross section
  if (fabs(Sy) + fabs(Sz) > 1.0e-12)
    shearTransform(EK,Sy,Sz);

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
                                         double Aly, double Alz, double ItoA,
                                         double Sy,  double Sz) const
{
  double L2 = L*L;
  double Cy = N/(L*(1.0+Aly)*(1.0+Aly));
  double Cz = N/(L*(1.0+Alz)*(1.0+Alz));

  EK.resize(12,12,true);
  EK( 2, 2) =  EK( 8, 8) =  Cy*(1.2+Aly*(2.0+Aly));
  EK( 3, 3) =  EK( 9, 9) =  Cz*(1.2+Alz*(2.0+Alz));
  EK( 4, 4) =  EK(10,10) =  N*ItoA/L;
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
  for (size_t i = 2; i <= 12; i++)
    for (size_t j = 1; j < i; j++)
      EK(i,j) = EK(j,i);

  // Adjustment due to non-symmetric cross section
  if (fabs(Sy) + fabs(Sz) > 1.0e-12)
    shearTransform(EK,Sy,Sz);

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
  for (size_t i = 2; i <= 12; i++)
    for (size_t j = 1; j < i; j++)
      EM(i,j) = EM(j,i);

#if INT_DEBUG > 1
  std::cout <<"ElasticBeam: local mass matrix:"<< EM;
#endif
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
  double EA, EIy, EIz, GIt, Aly, Alz, ItoA;
  double rhoA, I_xx, I_yy, I_zz, CG_y, CG_z, Sy, Sz;
  bool hasGrF = gravity.isZero() ? false : eS > 0;
  myProp->eval(X, L0, E, G, rho, hasGrF, eM > 0,
               EA, EIy, EIz, GIt, Aly, Alz,
               rhoA, CG_y, CG_z, I_xx, I_yy, I_zz, ItoA, Sy, Sz);
#if INT_DEBUG > 1
  std::cout <<"\n             EA = "<< EA
            <<" EI = "<< EIy <<" "<< EIz <<" GIt = "<< GIt
            <<", Alpha_y = "<< Aly <<" Alpha_z = "<< Alz;
  std::cout <<", rho*A = "<< rhoA <<" rho*I = "<< I_xx <<" "<< I_yy <<" "<< I_zz
            <<", CoG = "<< CG_y <<" "<< CG_z
            <<", Sy = "<< Sy <<" Sz = "<< Sz << std::endl;
#endif

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  Matrix tmpEK;

  if (eKm) // Evaluate the material stiffness matrix
    this->getMaterialStiffness(elMat.A[eKm-1],L0,EA,GIt,EIy,EIz,Aly,Alz,Sy,Sz);

  Vector v;
  double N = 0.0;

  if ((iS || eKg) && eV.normInf() > 1.0e-16*L0)
  {
    v = eV; // Transform the element displacement vector to local coordinates
    if (!utl::transform(v,this->getLocalAxes(elmInt),true))
      return false;

#if INT_DEBUG > 1
    std::cout <<"ElasticBeam: v"<< v;
#endif
    if (eKg) N = EA*(v(7)-v(1))/L0; // Axial force
  }

  bool hasIntS = v.empty() ? false : iS > 0;
  if (hasIntS)
  {
    v *= -1.0;

    // Calculate internal forces, S_int = Km*v
    if (!eKm) this->getMaterialStiffness(tmpEK,L0,EA,GIt,EIy,EIz,Aly,Alz,Sy,Sz);
    Matrix& Km = eKm ? elMat.A[eKm-1] : tmpEK;
    if (!Km.multiply(v,elMat.b[iS-1]))
      return false;

#if INT_DEBUG > 1
    std::cout <<"ElasticBeam: -S_int"<< elMat.b[iS-1];
#endif
  }

  if (eKg && N != 0.0)
  {
#if INT_DEBUG > 1
    std::cout <<"ElasticBeam: Axial force, N = "<< N;
#endif

    // Evaluate the geometric stiffness matrix
    Matrix& Kg = eKg == eKm ? tmpEK : elMat.A[eKg-1];
    this->getGeometricStiffness(Kg,N,L0,EIy,EIz,Aly,Alz,ItoA,Sy,Sz);
    if (eKg == eKm) elMat.A[eKm-1].add(Kg);
  }

  Matrix& Mm = eM ? elMat.A[eM-1] : tmpEK;
  if (eM || hasGrF)
  {
    // Evaluate the mass matrix
    this->getMassMatrix(Mm,rhoA,I_xx,I_yy,I_zz,L0);
    if (CG_y != 0.0 || CG_z != 0.0) // Transform to neutral axis location
      eccTransform(Mm,Vec3(0.0,-CG_y,-CG_z),Vec3(0.0,-CG_y,-CG_z));
  }

  if (hasGrF)
  {
    // Calculate external (gravitation) forces, S_ext = M*g
    Vec3 gvec = gravity*this->getLocalAxes(elmInt);
    Vector& S = elMat.b[eS-1];
    if (!hasIntS) S.fill(0.0);
    for (size_t i = 1; i <= S.size(); i++)
      for (size_t j = 0; j < Mm.cols(); j += npv)
        S(i) += Mm(i,j+1)*gvec.x + Mm(i,j+2)*gvec.y + Mm(i,j+3)*gvec.z;

#if INT_DEBUG > 1
    std::cout <<"ElasticBeam: "<< (eS == iS ? "S_ext - S_int" : "S_ext") << S;
#endif
  }

  return true;
}


/*!
  This method is overridden to also transform the element matrices from
  the local axes of the beam element to the global coordinate axes.
  This includes the effects of eccentric end points, if any.
*/

bool ElasticBeam::finalizeElement (LocalIntegral& elmInt,
                                   const TimeDomain& time, size_t)
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (inLocalAxes)
  {
    const Matrix& Tlg = this->getLocalAxes(elmInt);

    // Transform the element matrices to global coordinates
    for (Matrix& A : elMat.A)
      if (!utl::transform(A,Tlg))
        return false;

    // Transform the element force vectors to global coordinates
    for (Vector& b : elMat.b)
      if (!utl::transform(b,Tlg))
        return false;
  }

  const Vec3& ecc1 = myProp->ecc1; // End 1 offset in global coordinates
  const Vec3& ecc2 = myProp->ecc2; // End 2 offset in global coordinates

  if (!ecc1.isZero() || !ecc2.isZero())
  {
    // Transform all element matrices from eccentric end points to grid points
    for (Matrix& A : elMat.A)
      eccTransform(A,ecc1,ecc2);
    for (Vector& b : elMat.b)
      eccTransform(b,ecc1,ecc2);
  }

  return this->ElasticBase::finalizeElement(elmInt,time);
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
  return fld == 0 || !anasol ? 1 : 3;
}


std::string ElasticBeamNorm::getName (size_t, size_t j, const char* prfix) const
{
  static const char* u[3] = {
     "||u^h||_L2",
     "||u||_L2",
     "||e||_L2, e=u^h-u"
  };

  return prfix ? prfix + std::string(" ") + u[j-1] : u[j-1];
}


/*!
  This just forwards to the base class method not using any FiniteElement data,
  since the norm evaluation does not need the local-to-global transformation.
*/

bool ElasticBeamNorm::initElement (const std::vector<int>& MNPC,
                                   const FiniteElement&, const Vec3&, size_t,
                                   LocalIntegral& elmInt)
{
  return this->initElement(MNPC,elmInt);
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
