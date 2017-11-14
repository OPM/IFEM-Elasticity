// $Id$
//==============================================================================
//!
//! \file SIMLinEl3D.C
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for 3D NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#include "SIMLinEl.h"
#include "AnalyticSolutions.h"
#include "Vec3Oper.h"

//! Option for output of Gauss points to VTF for 3D problems.
template<> bool SIMLinEl3D::GIpointsVTF = false;


template<> bool SIMLinEl3D::parseDimSpecific (char* keyWord, std::istream& is)
{
  if (!strncasecmp(keyWord,"ANASOL",6)) {
    int code = -1;
    char* cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"HOLE",4))
    {
      double a  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      double nu = atof(strtok(NULL," "));
      IFEM::cout <<"\nAnalytical solution: Hole a="<< a <<" F0="<< F0
                 <<" nu="<< nu << std::endl;
      if (!mySol)
        mySol = new AnaSol(new Hole(a,F0,nu,true));
    }
    else if (!strncasecmp(cline,"LSHAPE",6))
    {
      double a  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      double nu = atof(strtok(NULL," "));
      IFEM::cout <<"\nAnalytical solution: Lshape a="<< a <<" F0="<< F0
                 <<" nu="<< nu << std::endl;
      if (!mySol)
        mySol = new AnaSol(new Lshape(a,F0,nu,true));
    }
    else if (!strncasecmp(cline,"CANTS",5))
    {
      double L  = atof(strtok(NULL," "));
      double H  = atof(strtok(NULL," "));
      double F0 = atof(strtok(NULL," "));
      IFEM::cout <<"\nAnalytical solution: CanTS L="<< L <<" H="<< H
                 <<" F0="<< F0 << std::endl;
      if (!mySol)
        mySol = new AnaSol(new CanTS(L,H,F0,true));
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      IFEM::cout <<"\nAnalytical solution: Expression"<< std::endl;
      int lines = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      if (!mySol)
        mySol = new AnaSol(is,lines,false);
    }
    else
    {
      std::cerr <<"  ** SIMLinEl3D::parse: Invalid analytical solution "
                << cline <<" (ignored)"<< std::endl;
      return true;
    }

    // Define the analytical boundary traction field
    if (code == -1)
      code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
    if (code > 0 && mySol->getStressSol())
    {
      IFEM::cout <<"Pressure code "<< code
                 <<": Analytical traction"<< std::endl;
      this->setPropertyType(code,Property::NEUMANN);
      myTracs[code] = new TractionField(*mySol->getStressSol());
    }
  }

  return false;
}


template<> bool SIMLinEl3D::parseDimSpecific (const TiXmlElement* child)
{
  if (!strcasecmp(child->Value(),"anasol")) {
    std::string type;
    utl::getAttribute(child,"type",type,true);
    if (type == "hole") {
      double a = 0.0, F0 = 0.0, nu = 0.0;
      utl::getAttribute(child,"a",a);
      utl::getAttribute(child,"F0",F0);
      utl::getAttribute(child,"nu",nu);
      IFEM::cout <<"\tAnalytical solution: Hole a="<< a <<" F0="<< F0
                 <<" nu="<< nu << std::endl;
      if (!mySol)
        mySol = new AnaSol(new Hole(a,F0,nu,true));
    }
    else if (type == "lshape") {
      double a = 0.0, F0 = 0.0, nu = 0.0;
      utl::getAttribute(child,"a",a);
      utl::getAttribute(child,"F0",F0);
      utl::getAttribute(child,"nu",nu);
      IFEM::cout <<"\tAnalytical solution: Lshape a="<< a <<" F0="<< F0
                 <<" nu="<< nu << std::endl;
      if (!mySol)
        mySol = new AnaSol(new Lshape(a,F0,nu,true));
    }
    else if (type == "cants") {
      double L = 0.0, H = 0.0, F0 = 0.0;
      utl::getAttribute(child,"L",L);
      utl::getAttribute(child,"H",H);
      utl::getAttribute(child,"F0",F0);
      IFEM::cout <<"\tAnalytical solution: CanTS L="<< L <<" H="<< H
                 <<" F0="<< F0 << std::endl;
      if (!mySol)
        mySol = new AnaSol(new CanTS(L,H,F0,true));
    }
    else if (type == "pipe") {
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
        mySol = new AnaSol(new Pipe(Ri,Ro,Ti,To,T0,E,nu,alpha,true,polar));
    }
    else if (type == "expression" || type == "fields") {
      type[0] = toupper(type[0]);
      IFEM::cout <<"\tAnalytical solution: "<< type << std::endl;
      if (!mySol)
        mySol = new AnaSol(child,false);
    }
    else
      std::cerr <<"  ** SIMLinEl3D::parse: Invalid analytical solution "
                << type <<" (ignored)"<< std::endl;

    // Define the analytical boundary traction field
    int code = 0;
    utl::getAttribute(child,"code",code);
    if (code > 0 && mySol && mySol->getStressSol()) {
      IFEM::cout <<"\tNeumann code "<< code
                 <<": Analytical traction"<< std::endl;
      setPropertyType(code,Property::NEUMANN);
      myTracs[code] = new TractionField(*mySol->getStressSol());
    }
  }
  else
    return false;

  return true;
}
