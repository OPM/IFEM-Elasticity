// $Id$
//==============================================================================
//!
//! \file LocalSystems.C
//!
//! \date Apr 07 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Local coordinate systems for linear elasticity problems.
//!
//==============================================================================

#include "IFEM.h"
#include "Elasticity.h"
#include "Utilities.h"
#include "Tensor.h"
#include "tinyxml.h"
#ifdef PRINT_CS
#include <fstream>
#endif


/*!
  \brief Local coordinate system for a cylinder along global z-axis.
*/

class CylinderCS : public LocalSystem
{
public:
  //! \brief The constructor prints a message making user aware of its presense.
  CylinderCS()
  {
    IFEM::cout <<"\nLocal coordinate system: Cylindric"<< std::endl;
  }

  //! \brief Empty destructor.
  virtual ~CylinderCS() {}

  //! \brief Computes the global-to-local transformation at the point \a X.
  virtual const Tensor& getTmat(const Vec3& X) const
  {
    static Tensor T(3);
    double r = hypot(X.x,X.y);
    T(1,1) = X.x/r;
    T(1,2) = X.y/r;
    T(2,1) = -T(1,2);
    T(2,2) = T(1,1);
    T(3,3) = 1.0;
    return T;
  }
};


/*!
  \brief Local coordinate system for a cylinder along the global z-axis,
  closed by a spherical cap.
*/

class CylinderSphereCS : public LocalSystem
{
public:
  //! \brief The constructor prints a message making user aware of its presense.
  CylinderSphereCS(double H = 0.0) : h(H)
  {
    IFEM::cout <<"\nLocal coordinate system: Cylindric with Spherical cap, h="
               << h << std::endl;
#ifdef PRINT_CS
    sn.open("nodes.dat");
    se.open("elements.dat");
    s1.open("v1.dat");
    s2.open("v2.dat");
    s3.open("v3.dat");
    sn <<"\n*NODES 4\n";
    se <<"\n*ELEMENTS 4\n%NODES #4\n"
       <<"%NO_ID\n%MAP_NODE_INDICES\n%PART_ID 4\n%POINTS\n";
    s1 <<"\n*RESULTS 31\n%NO_ID\n%DIMENSION 3\n%PER_NODE #4\n";
    s2 <<"\n*RESULTS 32\n%NO_ID\n%DIMENSION 3\n%PER_NODE #4\n";
    s3 <<"\n*RESULTS 33\n%NO_ID\n%DIMENSION 3\n%PER_NODE #4\n";
  }

  //! \brief The destructor finalizes the VTF files plotting the local axes.
  virtual ~CylinderSphereCS()
  {
    s1 <<"\n*GLVIEWVECTOR 2\n%NAME \"v1\"\n%STEP 1\n31\n";
    s2 <<"\n*GLVIEWVECTOR 3\n%NAME \"v2\"\n%STEP 1\n32\n";
    s3 <<"\n*GLVIEWVECTOR 4\n%NAME \"v3\"\n%STEP 1\n33\n";
#endif
  }

  //! \brief Computes the global-to-local transformation at the point \a X.
  virtual const Tensor& getTmat(const Vec3& X) const
  {
    static Tensor T(3);
#ifdef PRINT_CS
    sn << X <<'\n';
    static int iel = 0;
    se << ++iel <<'\n';
#endif
    if (patch == 1) // Cylindric system {-z,theta,r}
    {
      T.zero();
      double r = hypot(X.x,X.y);
      T(1,3) = -1.0;
      T(2,1) = -X.y/r;
      T(2,2) =  X.x/r;
      T(3,1) =  T(2,2);
      T(3,2) = -T(2,1);
#ifdef PRINT_CS
      s1 <<"0 0 -1\n";
      s2 << T(2,1) <<" "<< T(2,2) <<" 0\n";
      s3 << T(3,1) <<" "<< T(3,2) <<" 0\n";
#endif
    }
    else // Spherical system {phi,theta,r}
    {
      Vec3 v3(X.x,X.y,X.z-h);
      v3 /= v3.length();
      double theta = atan2(X.y,X.x);
      double phi = acos(v3.z);
      Vec3 v1(cos(theta)*cos(phi),sin(theta)*cos(phi),-sin(phi));
      Vec3 v2(v3,v1);
      for (int i = 1; i <= 3; i++)
      {
        T(1,i) = v1[i-1];
        T(2,i) = v2[i-1];
        T(3,i) = v3[i-1];
      }
#ifdef PRINT_CS
      s1 << v1 <<'\n';
      s2 << v2 <<'\n';
      s3 << v3 <<'\n';
#endif
    }
    return T;
  }

private:
  double h; //!< Height above global origin of the centre of the sphere
#ifdef PRINT_CS
  mutable std::ofstream sn; //!< VTF output stream for CS nodes
  mutable std::ofstream se; //!< VTF output stream for CS point elements
  mutable std::ofstream s1; //!< VTF output stream for vector v1 of local CS
  mutable std::ofstream s2; //!< VTF output stream for vector v2 of local CS
  mutable std::ofstream s3; //!< VTF output stream for vector v3 of local CS
#endif
};


bool Elasticity::parseLocalSystem (const char* cline)
{
  if (!strncasecmp(cline,"CYLINDRICZ",10))
    this->setLocalSystem(new CylinderCS);
  else if (!strncasecmp(cline,"CYLINDER+SPHERE",15))
    this->setLocalSystem(new CylinderSphereCS(atof(cline+15)));
  else
  {
    std::cerr <<"  ** Unsupported coordinate system: "
	      << cline <<" (ignored)"<< std::endl;
    return false;
  }

  return true;
}


bool Elasticity::parseLocalSystem (const TiXmlElement* elem)
{
  if (elem->FirstChild() == NULL) return false;

  // Caution: When running adaptively, the below will cause a small memory
  // leak because the coordinate system objects are only deleted by the
  // Elasticity destructor (and not in SIMbase::clearProperties).
  if (!strcasecmp(elem->FirstChild()->Value(),"cylindricz"))
    this->setLocalSystem(new CylinderCS);
  else if (!strcasecmp(elem->FirstChild()->Value(),"cylinder+sphere"))
  {
    double H = 0.0;
    utl::getAttribute(elem,"H",H);
    this->setLocalSystem(new CylinderSphereCS(H));
  }
  else
    std::cerr <<"  ** Unsupported local coordinate system: "
	      << elem->FirstChild()->Value() <<" (ignored)"<< std::endl;

  return true;
}
