// $Id$
//==============================================================================
//!
//! \file ElasticBeam.C
//!
//! \date Aug 10 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class representing a non-linear elastic beam.
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

//! Element matrices and data for non-linear dynamic beam FE problems
typedef BeamMats<HHTMats> HHTBeamMats;

template<> HHTBeamMats::BeamMats (double a1, double a2, double b, double)
  : HHTMats(a1,a2,b,true) {}


ElasticBeam::ElasticBeam (unsigned short int n) : inLocalAxes(true)
{
  npv = 6; // Number of primary unknowns per node
  nSV = n; // Number of solution vectors in core

  // Default material parameters
  E = 2.05e11;
  G = 7.94e10;
  rho = 7.85e3;

  // Default cross section parameters
  A  = 0.1;
  Ix = It = 0.002;
  Iy = Iz = 0.001;
  Ky = Kz = 0.0;
  Sy = Sz = 0.0;

  this->initPropFunc();
}


ElasticBeam::~ElasticBeam ()
{
  delete lineLoad;
  delete cplLoad;
  delete EAfunc;
  delete EIyfunc;
  delete EIzfunc;
  delete GItfunc;
  delete rhofunc;
  delete Ixfunc;
  delete Iyfunc;
  delete Izfunc;
  delete CGyfunc;
  delete CGzfunc;
}


void ElasticBeam::initPropFunc ()
{
  lineLoad = cplLoad = nullptr;
  EAfunc = EIyfunc = EIzfunc = GItfunc = nullptr;
  rhofunc = Ixfunc = Iyfunc = Izfunc = nullptr;
  CGyfunc = CGzfunc = nullptr;
}


void ElasticBeam::printLog () const
{
  IFEM::cout <<"ElasticBeam: E = "<< E <<", G = "<< G <<", rho = "<< rho
             <<"\n             A = "<< A <<" Ix = "<< Ix
             <<", Iy = "<< Iy <<", Iz = "<< Iz <<", It = "<< It
             <<"\n             Ky = "<< Ky <<", Kz = "<< Kz
             <<", Sy = "<< Sy <<", Sz = "<< Sz << std::endl;
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


bool ElasticBeam::parsePipe (const TiXmlElement* prop, double& A, double& I)
{
  double D, R;
  if (!utl::getAttribute(prop,"R",R))
  {
    if (utl::getAttribute(prop,"D",D))
      R = 0.5*D; // Radius of pipe cross section
    else
      return false;
  }

  double R2 = R*R, r2 = 0.0, t = 0.0;
  if (utl::getAttribute(prop,"t",t))
    r2 = (R-t)*(R-t); // Inner radius of hollow pipe

  A = M_PI*(R2-r2);
  I = M_PI*(R2*R2-r2*r2)*0.25;
  return true;
}


bool ElasticBeam::parseBox (const TiXmlElement* prop,
                            double& A, double& Iy, double& Iz)
{
  double H = 0.0;
  if (!utl::getAttribute(prop,"H",H))
    return false;

  double B = H;
  utl::getAttribute(prop,"B",B);

  A  = B*H;
  Iy = A*H*H/12.0;
  Iz = A*B*B/12.0;
  return true;
}


void ElasticBeam::parseBeamProperties (const TiXmlElement* prop)
{
  if (ElasticBeam::parsePipe(prop,A,Iz))
  {
    Iy = Iz;
    It = Ix = Iz*2.0;
    Ky = Kz = 2.0;
  }
  else if (ElasticBeam::parseBox(prop,A,Iy,Iz))
  {
    It = Ix = Iy + Iz;
    Ky = Kz = 1.2;
  }

  utl::getAttribute(prop,"A",A);
  utl::getAttribute(prop,"Ix",Ix);
  utl::getAttribute(prop,"Iy",Iy);
  utl::getAttribute(prop,"Iz",Iz);
  utl::getAttribute(prop,"It",It);
  utl::getAttribute(prop,"Sy",Sy);
  utl::getAttribute(prop,"Sz",Sz);
  utl::getAttribute(prop,"Ky",Ky);
  utl::getAttribute(prop,"Kz",Kz);
  IFEM::cout <<"    Constant beam properties:"
             <<"\n\tCross section area = "<< A
             <<", moments of inertia = "<< Ix <<" "<< Iy <<" "<< Iz <<" "<< It
             <<"\n\tShear parameters = "<< Sy <<" "<< Sz <<" "<< Ky <<" "<< Kz
             << std::endl;

  RealFunc** pf = nullptr;
  const TiXmlElement* child = prop->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (child->FirstChild())
    {
      if (!pf)
        IFEM::cout <<"    Continuous beam properties:\n";
      IFEM::cout <<"\t"<< child->Value();
      if (!strcmp(child->Value(),"EA"))
        pf = &EAfunc;
      else if (!strcmp(child->Value(),"EIy"))
        pf = &EIyfunc;
      else if (!strcmp(child->Value(),"EIz"))
        pf = &EIzfunc;
      else if (!strcmp(child->Value(),"GIt"))
        pf = &GItfunc;
      else if (!strcmp(child->Value(),"rho"))
        pf = &rhofunc;
      else if (!strcmp(child->Value(),"Ix"))
        pf = &Ixfunc;
      else if (!strcmp(child->Value(),"Iy"))
        pf = &Iyfunc;
      else if (!strcmp(child->Value(),"Iz"))
        pf = &Izfunc;
      else if (!strcmp(child->Value(),"CGy"))
        pf = &CGyfunc;
      else if (!strcmp(child->Value(),"CGz"))
        pf = &CGzfunc;
      else
      {
        IFEM::cout <<" (ignored)"<< std::endl;
        continue;
      }

      std::string type;
      utl::getAttribute(child,"type",type);
      if (!type.empty()) IFEM::cout <<" ("<< type <<")";
      *pf = utl::parseRealFunc(child->FirstChild()->Value(),type);
      IFEM::cout << std::endl;
    }
}


LocalIntegral* ElasticBeam::getLocalIntegral (size_t, size_t, bool) const
{
  ElmMats* result;
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
  if (e3*e1 < e2*e1)
  {
    e2.cross(e3,e1); // Local Y-axis = e3xe1 / |e3xe1|
    e2.normalize();
    e3.cross(e1,e2); // Local Z-axis = e1xe2
  }
  else
  {
    e3.cross(e1,e2); // Local Z-axis = e1xe2 / |e1xe2|
    e3.normalize();
    e2.cross(e3,e1); // Local Y-axis = e3xe1
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

void ElasticBeam::getMaterialStiffness (Matrix& EK, double EA, double GIt,
                                        double EIy, double EIz, double L) const
{
  double L2  = L*L;
  double L3  = L*L2;
  double Aly = 12.0*EIy*Ky/(G*A*L2);
  double Alz = 12.0*EIz*Kz/(G*A*L2);
#if INT_DEBUG > 1
  std::cout <<"ElasticBeam: Alpha_y = "<< Aly <<" Alpha_z = "<< Alz;
#endif

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

  if (fabs(Sy) + fabs(Sz) > 0.00001*sqrt(A))
  {
    // Adjustment due to non-symmetric cross section
    for (j = 1; j <= 12; j++)
      for (i = 0; i <= 6; i += 6)
        EK(i+4,j) -= Sz*EK(i+2,j) + Sy*EK(i+3,j);

    for (i = 1; i <= 12; i++)
      for (j = 0; j <= 6; j += 6)
        EK(i,j+4) -= EK(i,j+2)*Sz + EK(i,j+3)*Sy;
  }

#if INT_DEBUG > 1
  std::cout <<"\nElasticBeam: local material stiffness matrix:"<< EK;
#endif
}


/*!
  \details This matrix is taken from page 14 in
  http://people.duke.edu/~hpgavin/cee421/frame-finite-def.pdf
*/

void ElasticBeam::getGeometricStiffness (Matrix& EK, double EIy, double EIz,
                                         double L, double N) const
{
  double L2  = L*L;
  double Aly = 12.0*EIy*Ky/(G*A*L2);
  double Alz = 12.0*EIz*Kz/(G*A*L2);
  double Cy  = N/(L*(1.0+Aly)*(1.0+Aly));
  double Cz  = N/(L*(1.0+Alz)*(1.0+Alz));

  EK.resize(12,12,true);
  EK( 2, 2) =  EK( 8, 8) =  Cy*(1.2+Aly*(2.0+Aly));
  EK( 3, 3) =  EK( 9, 9) =  Cz*(1.2+Alz*(2.0+Alz));
  EK( 4, 4) =  EK(10,10) =  N*It/(A*L);
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

  if (fabs(Sy) + fabs(Sz) > 0.00001*sqrt(A))
  {
    // Adjustment due to non-symmetric cross section
    for (j = 1; j <= 12; j++)
      for (i = 0; i <= 6; i += 6)
        EK(i+4,j) -= Sz*EK(i+2,j) + Sy*EK(i+3,j);

    for (i = 1; i <= 12; i++)
      for (j = 0; j <= 6; j += 6)
        EK(i,j+4) -= EK(i,j+2)*Sz + EK(i,j+3)*Sy;
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
  double L, L0 = X0.length();
  if (L0 <= 1.0e-8)
  {
    std::cerr <<" *** ElasticBeam::evalInt: Zero initial element length "
              << L0 << std::endl;
    return false;
  }
#if INT_DEBUG > 1
  std::cout <<"ElasticBeam: X = "<< X <<", L0 = "<< L0;
#endif

  const Vector& eV = elmInt.vec.front();
  if (eV.empty())
    L = L0;
  else
  {
    // Calculate current element length
    Vec3 U1(eV.ptr(),npv);
    Vec3 U2(eV.ptr()+eV.size()-npv,npv);
    Vec3 X1 = X0 + U2 - U1;
    L = X1.length();
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

  // Evaluate beam stiffness properties at this point
  double EA  = EAfunc  ? (*EAfunc)(X)  : E*A;
  double EIy = EIyfunc ? (*EIyfunc)(X) : E*Iy;
  double EIz = EIzfunc ? (*EIzfunc)(X) : E*Iz;
  double GIt = GItfunc ? (*GItfunc)(X) : G*It;
#if INT_DEBUG > 1
  std::cout <<"\n             EA = "<< EA
            <<" EI = "<< EIy <<" "<< EIz <<" GIt = "<< GIt;
#endif

  // Evaluate the beam mass properties (if needed) at this point
  bool hasGrF = gravity.isZero() ? false : eS > 0;
  double rhoA = rhofunc && (eM > 0 || hasGrF) ? (*rhofunc)(X) : rho*A;
  double I_xx = Ixfunc  &&  eM > 0            ? (*Ixfunc)(X)  : rho*Ix;
  double I_yy = Iyfunc  &&  eM > 0            ? (*Iyfunc)(X)  : rho*Iy;
  double I_zz = Izfunc  &&  eM > 0            ? (*Izfunc)(X)  : rho*Iz;
  double CG_y = CGyfunc && (eM > 0 || hasGrF) ? (*CGyfunc)(X) : 0.0;
  double CG_z = CGzfunc && (eM > 0 || hasGrF) ? (*CGzfunc)(X) : 0.0;
#if INT_DEBUG > 1
  std::cout <<", rho*A = "<< rhoA <<" rho*I = "<< I_xx <<" "<< I_yy <<" "<< I_zz
            <<", CoG = "<< CG_y <<" "<< CG_z << std::endl;
#endif

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (hasGrF)
  {
    // External (gravitation) forces
    Vector& S = elMat.b[eS-1];
    Vec3 gvec = (0.5*rhoA*L0)*gravity; // Nodal gravity force at each end
    for (unsigned short int i = 1; i <= 3; i++)
      S(i) = S(npv+i) = gvec[i-1];

    if (CG_y != 0.0 || CG_z != 0.0)
    {
      // The centre of gravity has an offset w.r.t. the neutral axis and will
      // therefore result in an additional torque load on the element
      Tensor Tlg(this->getLocalAxes(elmInt)); // Local-to-global transformation
      Vec3   gloc = gvec * Tlg; // Nodal gravity force in local element axes
      double Tgrav = gloc.z*CG_y - gloc.y*CG_z; // Torque from eccentric gravity
      Vec3   gmom = Tgrav*Tlg[0]; // Global moment from the eccentric gravity
      for (unsigned short int i = 4; i <= 6; i++)
        S(i) = S(npv+i) = gmom[i-4];
    }
#if INT_DEBUG > 1
    std::cout <<"ElasticBeam: S_ext"<< S << std::endl;
#endif
  }

  if (eKm) // Evaluate the material stiffness matrix
    this->getMaterialStiffness(elMat.A[eKm-1],EA,GIt,EIy,EIz,L0);

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
    if (!eKm) this->getMaterialStiffness(tmpKm,EA,GIt,EIy,EIz,L0);
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
      this->getGeometricStiffness(Kg,EIy,EIz,L0,N);
      elMat.A[eKm-1].add(Kg);
    }
    else
      this->getGeometricStiffness(elMat.A[eKg-1],EIy,EIz,L0,N);
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
  ElasticBeam& problem = static_cast<ElasticBeam&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the FE displacements
  Vec3 Uh = problem.displacement(fe,elmInt.vec.front());

  // L2-norm of FE displacements
  size_t ip = 0;
  pnorm[ip++] = Uh*Uh*fe.detJxW;

  if (anasol)
  {
    // Exact displacements
    Vec3 U = (*anasol)(X);

    // L2-norm of exact solution
    pnorm[ip++] += U*U*fe.detJxW;

    // L2-norm of displacement error
    U -= Uh;
    pnorm[ip++] += U*U*fe.detJxW;
  }

  return true;
}
