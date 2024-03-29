// $Id$
//==============================================================================
//!
//! \file SIMLinEl2D.C
//!
//! \date Dec 04 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for 2D NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#include "SIMLinEl.h"
#include "AnalyticSolutions.h"


template<> bool SIMLinEl2D::parseDimSpecific (char* cline)
{
  if (!strncasecmp(cline,"HOLE",4))
  {
    double a  = atof(strtok(nullptr," "));
    double F0 = atof(strtok(nullptr," "));
    double nu = atof(strtok(nullptr," "));
    IFEM::cout <<"\nAnalytical solution: Hole a="<< a <<" F0="<< F0
               <<" nu="<< nu << std::endl;
    if (!mySol)
      mySol = new AnaSol(new Hole(a,F0,nu));
  }
  else if (!strncasecmp(cline,"LSHAPE",6))
  {
    double a  = atof(strtok(nullptr," "));
    double F0 = atof(strtok(nullptr," "));
    double nu = atof(strtok(nullptr," "));
    IFEM::cout <<"\nAnalytical solution: Lshape a="<< a <<" F0="<< F0
               <<" nu="<< nu << std::endl;
    if (!mySol)
      mySol = new AnaSol(new Lshape(a,F0,nu));
  }
  else if (!strncasecmp(cline,"CANTS",5))
  {
    double L  = atof(strtok(nullptr," "));
    double H  = atof(strtok(nullptr," "));
    double F0 = atof(strtok(nullptr," "));
    IFEM::cout <<"\nAnalytical solution: CanTS L="<< L <<" H="<< H
               <<" F0="<< F0 << std::endl;
    if (!mySol)
      mySol = new AnaSol(new CanTS(L,H,F0));
  }
  else if (!strncasecmp(cline,"CANTM",5))
  {
    double H  = atof(strtok(nullptr," "));
    double M0 = atof(strtok(nullptr," "));
    IFEM::cout <<"\nAnalytical solution: CanTM H="<< H
               <<" M0="<< M0 << std::endl;
    if (!mySol)
      mySol = new AnaSol(new CanTM(H,M0));
  }
  else if (!strncasecmp(cline,"CURVED",6))
  {
    double a  = atof(strtok(nullptr," "));
    double b  = atof(strtok(nullptr," "));
    double u0 = atof(strtok(nullptr," "));
    double E  = atof(strtok(nullptr," "));
    IFEM::cout <<"\nAnalytical solution: Curved Beam a="<< a <<" b="<< b
               <<" u0="<< u0 <<" E="<< E << std::endl;
    if (!mySol)
      mySol = new AnaSol(new CurvedBeam(u0,a,b,E));
  }
  else
    return false;

  return true;
}


template<> bool SIMLinEl2D::parseDimSpecific (const tinyxml2::XMLElement* child,
                                              const std::string& type)
{
  if (type == "hole")
  {
    double a = 0.0, F0 = 0.0, nu = 0.0;
    utl::getAttribute(child,"a",a);
    utl::getAttribute(child,"F0",F0);
    utl::getAttribute(child,"nu",nu);
    IFEM::cout <<"\tAnalytical solution: Hole a="<< a <<" F0="<< F0
               <<" nu="<< nu << std::endl;
    if (!mySol)
      mySol = new AnaSol(new Hole(a,F0,nu));
  }
  else if (type == "lshape")
  {
    double a = 0.0, F0 = 0.0, nu = 0.0;
    utl::getAttribute(child,"a",a);
    utl::getAttribute(child,"F0",F0);
    utl::getAttribute(child,"nu",nu);
    IFEM::cout <<"\tAnalytical solution: Lshape a="<< a <<" F0="<< F0
               <<" nu="<< nu << std::endl;
    if (!mySol)
      mySol = new AnaSol(new Lshape(a,F0,nu));
  }
  else if (type == "cants")
  {
    double L = 0.0, H = 0.0, F0 = 0.0;
    utl::getAttribute(child,"L",L);
    utl::getAttribute(child,"H",H);
    utl::getAttribute(child,"F0",F0);
    IFEM::cout <<"\tAnalytical solution: CanTS L="<< L <<" H="<< H
               <<" F0="<< F0 << std::endl;
    if (!mySol)
      mySol = new AnaSol(new CanTS(L,H,F0));
  }
  else if (type == "cantm")
  {
    double L = 0.0, H = 0.0, M0 = 0.0;
    utl::getAttribute(child,"L",L);
    utl::getAttribute(child,"H",H);
    utl::getAttribute(child,"M0",M0);
    IFEM::cout <<"\tAnalytical solution: CanTM H="<< H
               <<" M0="<< M0 << std::endl;
    if (!mySol)
      mySol = new AnaSol(new CanTM(H,M0));
  }
  else if (type == "curved")
  {
    double a = 0.0, b = 0.0, u0 = 0.0, E = 0.0;
    utl::getAttribute(child,"a",a);
    utl::getAttribute(child,"b",b);
    utl::getAttribute(child,"u0",u0);
    utl::getAttribute(child,"E",E);
    IFEM::cout <<"\tAnalytical solution: Curved Beam a="<< a <<" b="<< b
               <<" u0="<< u0 <<" E="<< E << std::endl;
    if (!mySol)
      mySol = new AnaSol(new CurvedBeam(u0,a,b,E));
  }
  else if (type == "pipe")
  {
    double Ri = 0.0, Ro = 0.0, Ti = 0.0, To = 0.0, T0 = 273.0;
    double E = 0.0, nu = 0.0, alpha = 0.0;
    bool polar = false;
    utl::getAttribute(child,"Ri",Ri);
    utl::getAttribute(child,"Ro",Ro);
    utl::getAttribute(child,"Ti",Ti);
    utl::getAttribute(child,"To",To);
    utl::getAttribute(child,"Tref",T0);
    utl::getAttribute(child,"E",E);
    utl::getAttribute(child,"nu",nu);
    utl::getAttribute(child,"alpha",alpha);
    utl::getAttribute(child,"polar",polar);
    IFEM::cout <<"\tAnalytical solution: Pipe Ri="<< Ri <<" Ro="<< Ro
               <<" Ti="<< Ti <<" To="<< To << std::endl;
    if (!mySol)
      mySol = new AnaSol(new Pipe(Ri,Ro,Ti,To,T0,E,nu,alpha,false,polar));
  }
  else
    return false;

  return true;
}
