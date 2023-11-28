// $Id$
//==============================================================================
//!
//! \file BeamProperty.C
//!
//! \date Apr 7 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class representing a beam cross section property.
//!
//==============================================================================

#include "BeamProperty.h"
#include "Functions.h"
#include "Utilities.h"
#include "Vec3.h"
#include "IFEM.h"
#include "tinyxml2.h"


BeamProperty::BeamProperty (const tinyxml2::XMLElement* prop)
{
  // Default cross section parameters
  A  = 0.1;
  Ix = It = 0.002;
  Iy = Iz = 0.001;
  Ky = Kz = 0.0;
  Sy = Sz = 0.0;

  EAfunc = EIyfunc = EIzfunc = GItfunc = nullptr;
  rhofunc = Ixfunc = Iyfunc = Izfunc = nullptr;
  CGyfunc = CGzfunc = nullptr;

  if (prop)
    this->parse(prop);
}


BeamProperty::~BeamProperty ()
{
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


bool BeamProperty::parsePipe (const tinyxml2::XMLElement* prop, double& A, double& I)
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


bool BeamProperty::parseBox (const tinyxml2::XMLElement* prop,
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


void BeamProperty::parse (const tinyxml2::XMLElement* prop)
{
  if (parsePipe(prop,A,Iz))
  {
    Iy = Iz;
    It = Ix = Iz*2.0;
    Ky = Kz = 2.0;
  }
  else if (parseBox(prop,A,Iy,Iz))
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
  const tinyxml2::XMLElement* child = prop->FirstChildElement();
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


void BeamProperty::eval (const Vec3& X, double L,
                         double E, double G, double rho,
                         bool hasGrav, bool hasMass,
                         double& EA,   double& EI_y, double& EI_z,
                         double& GI_t, double& Al_y, double& Al_z,
                         double& rhoA, double& CG_y, double& CG_z,
                         double& I_xx, double& I_yy, double& I_zz) const
{
  // Evaluate beam stiffness properties at this point
  EA   = EAfunc  ? (*EAfunc)(X)  : E*A;
  EI_y = EIyfunc ? (*EIyfunc)(X) : E*Iy;
  EI_z = EIzfunc ? (*EIzfunc)(X) : E*Iz;
  GI_t = GItfunc ? (*GItfunc)(X) : G*It;
  Al_y = 12.0*EI_y*Ky/(G*A*L*L);
  Al_z = 12.0*EI_z*Kz/(G*A*L*L);

  // Evaluate the beam mass properties (if needed) at this point
  rhoA = rhofunc && (hasMass || hasGrav) ? (*rhofunc)(X) : rho*A;
  CG_y = CGyfunc && (hasMass || hasGrav) ? (*CGyfunc)(X) : 0.0;
  CG_z = CGzfunc && (hasMass || hasGrav) ? (*CGzfunc)(X) : 0.0;
  I_xx = Ixfunc  &&  hasMass             ? (*Ixfunc)(X)  : rho*Ix;
  I_yy = Iyfunc  &&  hasMass             ? (*Iyfunc)(X)  : rho*Iy;
  I_zz = Izfunc  &&  hasMass             ? (*Izfunc)(X)  : rho*Iz;
}


void BeamProperty::eval (const Vec3& X,
                         double E, double G, double rho,
                         double& EA, double& EI_y, double& EI_z,
                         double& GI_t,
                         double& rhoA, double& I_yy, double& I_zz) const
{
  // Evaluate beam stiffness properties at this point
  EA   = EAfunc  ? (*EAfunc)(X)  : E*A;
  EI_y = EIyfunc ? (*EIyfunc)(X) : E*Iy;
  EI_z = EIzfunc ? (*EIzfunc)(X) : E*Iz;
  GI_t = GItfunc ? (*GItfunc)(X) : G*It;

  // Evaluate the beam mass properties at this point
  rhoA = rhofunc ? (*rhofunc)(X) : rho*A;
  I_yy = Iyfunc  ? (*Iyfunc)(X)  : rho*Iy;
  I_zz = Izfunc  ? (*Izfunc)(X)  : rho*Iz;
}


double BeamProperty::evalRho (const Vec3& X, double rho) const
{
  return rhofunc ? (*rhofunc)(X) : rho*A;
}
