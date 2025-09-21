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
#include "Vec3Oper.h"
#include "IFEM.h"
#include "tinyxml2.h"
#include <fstream>
#include <cstring>


BeamProperty::BeamProperty (const tinyxml2::XMLElement* prop)
{
  // Default cross section parameters
  A  = 0.1;
  Ix = It = 0.002;
  Iy = Iz = 0.001;
  Ky = Kz = 0.0;
  Sy = Sz = 0.0;
  phi = 0.0;

  EAfunc = EIyfunc = EIzfunc = GItfunc = nullptr;
  rhofunc = Ixfunc = Iyfunc = Izfunc = nullptr;
  Axfunc = Ayfunc = Azfunc = nullptr;
  CGyfunc = CGzfunc = nullptr;
  Syfunc = Szfunc = nullptr;

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
  delete Axfunc;
  delete Ayfunc;
  delete Azfunc;
  delete Ixfunc;
  delete Iyfunc;
  delete Izfunc;
  delete CGyfunc;
  delete CGzfunc;
  delete Syfunc;
  delete Szfunc;
}


void BeamProperty::setConstant (const std::vector<double>& values)
{
  if (values.size() > 0) A   = values[0];
  if (values.size() > 2) Ix  = values[1] + values[2];
  if (values.size() > 1) Iy  = values[1];
  if (values.size() > 2) Iz  = values[2];
  if (values.size() > 3) It  = values[3];
  if (values.size() > 4) Ky  = values[4];
  if (values.size() > 5) Kz  = values[5];
  if (values.size() > 6) Sy  = values[6];
  if (values.size() > 7) Sz  = values[7];
  if (values.size() > 8) phi = values[8];
}


bool BeamProperty::parsePipe (const tinyxml2::XMLElement* prop,
                              double& A, double& I, double* K)
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

  if (K) *K = r2 > 0.0 ? 2.0 : 1.0/0.9;

  IFEM::cout <<"\tPipe("<< R <<","<< t
             <<"): A = "<< A <<" I = "<< I << std::endl;
  return true;
}


bool BeamProperty::parseBox (const tinyxml2::XMLElement* prop,
                             double& A, double& I1, double& I2,
                             double& J, double* K1, double* K2)
{
  double H = 0.0;
  if (!utl::getAttribute(prop,"H",H))
    return false;

  double B = H;
  utl::getAttribute(prop,"B",B);
  IFEM::cout <<"\tBox("<< H <<","<< B;

  double t1 = 0.0;
  if (utl::getAttribute(prop,"t",t1) || utl::getAttribute(prop,"t1",t1))
  {
    double t2 = t1;
    IFEM::cout <<","<< t1;
    if (utl::getAttribute(prop,"t2",t2))
      IFEM::cout <<","<< t2;

    double Bi = B - 2.0*t2;
    double Hi = H - 2.0*t1;

    A  = B*H - Bi*Hi;
    I1 = (B*H*H*H - Bi*Hi*Hi*Hi)/12.0;
    I2 = (H*B*B*B - Hi*Bi*Bi*Bi)/12.0;
    J  = 2.0*t2*t1*pow((B-t2)*(H-t1),2.0) / (B*t2 + H*t1 - t2*t2 - t1*t1);
    if (K1) *K1 = 0.5*A/(Bi*t1);
    if (K2) *K2 = 0.5*A/(Hi*t2);
  }
  else
  {
    A  = B*H;
    I1 = B*H*H*H/12.0;
    I2 = H*B*B*B/12.0;
    if (H > B) std::swap(B,H);
    J  = B*H*H*H*(1.0 - 0.63*(H/B)*(1.0-pow(H/B,4.0)/12.0))/3.0;
    if (K1) *K1 = 1.2;
    if (K2) *K2 = 1.2;
  }

  IFEM::cout <<"): A = "<< A <<" I = "<< I1 <<" "<< I2 <<" "<< J;
  if (K1 && K2) IFEM::cout <<" K = "<< *K1 <<" "<< *K2;
  IFEM::cout << std::endl;
  return true;
}


void BeamProperty::parse (const tinyxml2::XMLElement* prop)
{
  if (parsePipe(prop,A,Iz,&Kz))
  {
    Ky = Kz;
    Iy = Iz;
    It = Ix = Iz*2.0;
  }
  else if (parseBox(prop,A,Iy,Iz,It,&Ky,&Kz))
    Ix = Iy + Iz;

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

  if (prop->FirstChild() && this->readCSV(prop->FirstChild()->Value()))
    return; // Property functions are defined from a csv data file

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


bool BeamProperty::readCSV (const char* fileName)
{
  struct DataLine
  {
    std::string descr;
    std::string label;
    std::string unit;
    std::vector<double> values;
  } dline;

  std::ifstream fs(fileName);
  if (!fs) return false;

  std::map<std::string,DataLine> crossSections;
  char currentLine[BUFSIZ]; currentLine[0] = '#';
  while (fs.getline(currentLine+1,BUFSIZ-1))
  {
    dline.values.clear();
    char* val   = strtok(currentLine,",");
    dline.descr =  val ? val+1 : "";
    dline.label = (val = strtok(NULL,",")) ? val : "";
    dline.unit  = (val = strtok(NULL,",")) ? val : "";
    while ((val = strtok(NULL,",")))
    {
      dline.values.push_back(atof(val));
      // Scale to SI units [m]
      if (dline.unit == "[cm]")
        dline.values.back() *= 1.0e-2;
      else if (dline.unit == "[mm]")
        dline.values.back() *= 1.0e-3;
      else if (dline.unit == "[cm2]")
        dline.values.back() *= 1.0e-4;
      else if (dline.unit == "[mm2]")
        dline.values.back() *= 1.0e-6;
      else if (dline.unit == "[cm3]")
        dline.values.back() *= 1.0e-6;
      else if (dline.unit == "[mm3]")
        dline.values.back() *= 1.0e-9;
      else if (dline.unit == "[cm4]")
        dline.values.back() *= 1.0e-8;
      else if (dline.unit == "[mm4]")
        dline.values.back() *= 1.0e-12;
    }
    if (dline.unit.substr(1,2) == "cm" || dline.unit.substr(1,2) == "mm")
      dline.unit.erase(1,1);
#ifdef INT_DEBUG
    if (crossSections.empty())
      std::cout <<"Data from CSV-file "<< fileName <<":"<< std::endl;
    std::cout << dline.descr <<"\t"<< dline.label <<"\t"<< dline.unit;
    for (double dv : dline.values) std::cout <<" "<< dv;
    std::cout << std::endl;
#endif
    if (crossSections.find(dline.label) == crossSections.end())
      crossSections[dline.label] = dline;
  }

  if (crossSections.empty()) return false;

  // Lambda function returning the cross section data values for given key.
  auto&& getValues = [&crossSections,fileName](const std::string& key)
  {
    std::map<std::string,DataLine>::const_iterator it = crossSections.find(key);
    if (it != crossSections.end()) return it->second.values;

    std::cerr <<"  ** BeamProperty::readCSV: Property \""<< key
              <<"\" not found in "<< fileName << std::endl;
    static const std::vector<double> empty;
    return empty;
  };

  const std::vector<double>& xVal = getValues("x");
  if (xVal.empty()) return false; // Probably an invalid CSV file

#ifndef INT_DEBUG
  IFEM::cout <<"    Beam properties from CSV file "<< fileName <<":";
  for (const std::pair<const std::string,DataLine>& cs : crossSections)
    IFEM::cout <<"\n\t"<< cs.second.label <<" "<< cs.second.unit
               <<" "<< cs.second.descr;
  IFEM::cout << std::endl;
#endif

  // Lambda function creating a piece-wise linear property function.
  auto&& propertyFunc = [&xVal,getValues](const std::string& key)
  {
    const std::vector<double>& fVals = getValues(key);
    return fVals.empty() ? nullptr : new Interpolate1D(xVal,fVals);
  };

  Axfunc = propertyFunc("Ax");
  Ayfunc = propertyFunc("Ay");
  Azfunc = propertyFunc("Az");

  Ixfunc = propertyFunc("Ix");
  Iyfunc = propertyFunc("Iy");
  Izfunc = propertyFunc("Iz");

  Syfunc = propertyFunc("ey");
  Szfunc = propertyFunc("ez");

  return true;
}


void BeamProperty::eval (const Vec3& X, double L,
                         double E, double G, double rho,
                         bool hasGrav, bool hasMass,
                         double& EA,   double& EI_y, double& EI_z,
                         double& GI_t, double& Al_y, double& Al_z,
                         double& rhoA, double& CG_y, double& CG_z,
                         double& I_xx, double& I_yy, double& I_zz,
                         double& ItoA, double& S_y,  double& S_z) const
{
  // Evaluate beam stiffness properties at this point
  double Area = Axfunc ? (*Axfunc)(X) : A;
  EA   = EAfunc  ? (*EAfunc)(X)  : E*Area;
  EI_y = EIyfunc ? (*EIyfunc)(X) : E*(Iyfunc ? (*Iyfunc)(X) : Iy);
  EI_z = EIzfunc ? (*EIzfunc)(X) : E*(Izfunc ? (*Izfunc)(X) : Iz);
  GI_t = GItfunc ? (*GItfunc)(X) : G*(Ixfunc ? (*Ixfunc)(X) : It);
  Al_y = 12.0*(EI_y/(G*L*L)) * (Ayfunc ? 1.0 / (*Ayfunc)(X) : Ky / Area);
  Al_z = 12.0*(EI_z/(G*L*L)) * (Azfunc ? 1.0 / (*Azfunc)(X) : Kz / Area);
  ItoA = (Ixfunc ? (*Ixfunc)(X) : It) / Area;

  S_y = Syfunc ? (*Syfunc)(X) : Sy;
  S_z = Szfunc ? (*Szfunc)(X) : Sz;
  if (S_y*S_y + S_z*S_z < 1.0e-8*Area)
    S_y = S_z = 0.0;

  // Evaluate the beam mass properties (if needed) at this point
  rhoA = rhofunc && (hasMass || hasGrav) ? (*rhofunc)(X) : rho*Area;
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
  EA   = EAfunc  ? (*EAfunc)(X)  : E*(Axfunc ? (*Axfunc)(X) : A);
  EI_y = EIyfunc ? (*EIyfunc)(X) : E*(Iyfunc ? (*Iyfunc)(X) : Iy);
  EI_z = EIzfunc ? (*EIzfunc)(X) : E*(Izfunc ? (*Izfunc)(X) : Iz);
  GI_t = GItfunc ? (*GItfunc)(X) : G*(Ixfunc ? (*Ixfunc)(X) : It);

  // Evaluate the beam mass properties at this point
  rhoA = rhofunc ? (*rhofunc)(X) : rho*(Axfunc ? (*Axfunc)(X) : A);
  I_yy = Iyfunc  ? (*Iyfunc)(X)  : rho*(Iyfunc ? (*Iyfunc)(X) : Iy);
  I_zz = Izfunc  ? (*Izfunc)(X)  : rho*(Izfunc ? (*Izfunc)(X) : Iz);
}


double BeamProperty::evalRho (const Vec3& X, double rho) const
{
  return rhofunc ? (*rhofunc)(X) : rho*(Axfunc ? (*Axfunc)(X) : A);
}


std::ostream& operator<< (std::ostream& os, const BeamProperty& prop)
{
  const char* newline = "\n             ";

  os <<", A = "<< prop.A
     << newline <<"Ix = "<< prop.Ix <<", Iy = "<< prop.Iy <<", Iz = "<< prop.Iz
     <<", It = "<< prop.It << newline <<"Ky = "<< prop.Ky <<", Kz = "<< prop.Kz
     <<", Sy = "<< prop.Sy <<", Sz = "<< prop.Sz;
  if (fabs(prop.phi) > 1.0e-6) os << newline <<"phi = "<< prop.phi;
  if (!prop.ecc1.isZero()) os << newline <<"ecc1 = "<< prop.ecc1;
  if (!prop.ecc2.isZero()) os << newline <<"ecc2 = "<< prop.ecc2;

  return os;
}
