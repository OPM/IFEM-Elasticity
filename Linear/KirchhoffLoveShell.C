// $Id$
//==============================================================================
//!
//! \file KirchhoffLoveShell.C
//!
//! \date Feb 25 2018
//!
//! \author ... and ... / NTNU
//!
//! \brief Class for linear Kirchhoff-Love thin shell problems.
//!
//==============================================================================

#include "KirchhoffLoveShell.h"
#include "LinIsotropic.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "Utilities.h"
#include "IFEM.h"


void KirchhoffLoveShell::printLog () const
{
  IFEM::cout <<"KirchhoffLoveShell: thickness = "<< thickness
             <<", gravity = "<< gravity << std::endl;

  if (!material)
  {
    static LinIsotropic defaultMat;
    const_cast<KirchhoffLoveShell*>(this)->material = &defaultMat;
  }

  material->printLog();
}


bool KirchhoffLoveShell::evalInt (LocalIntegral& elmInt,
                                  const FiniteElement& fe,
                                  const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (eM) // Integrate the mass matrix
    ; // TODO: Include the mass-matrix terms in elmMat.A[eM-1]

  if (eK) // Integrate the stiffness matrix
    ; // TODO: Include the stiffness-matrix terms in elmMat.A[eK-1]

  if (eS) // Integrate the load vector due to gravitation and other body forces
    ; // TODO: Include the pressure/gravity load terms in elmMat.b[eS-1]

  return true;
}


bool KirchhoffLoveShell::evalBou (LocalIntegral& elmInt,
                                  const FiniteElement& fe,
                                  const Vec3& X, const Vec3& normal) const
{
  // TODO (if you want to support Neumann boundary conditions)
  std::cerr <<" *** KirchhoffLoveShell::evalBou not implemented."<< std::endl;
  return false;
}


bool KirchhoffLoveShell::evalSol (Vector& s,
                                  const FiniteElement& fe, const Vec3& X,
                                  const std::vector<int>& MNPC) const
{
  // Extract element displacements
  Vector eV;
  if (!primsol.empty() && !primsol.front().empty())
  {
    int ierr = utl::gather(MNPC,3,primsol.front(),eV);
    if (ierr > 0)
    {
      std::cerr <<" *** KirchhoffLoveShell::evalSol: Detected "
                << ierr <<" node numbers out of range."<< std::endl;
      return false;
    }
  }

  // Evaluate the stress resultant tensor
  return this->evalSol(s,eV,fe,X,true);
}


bool KirchhoffLoveShell::evalSol (Vector& s, const Vector& eV,
                                  const FiniteElement& fe, const Vec3& X,
                                  bool toLocal) const
{
  // TODO (if you want to support postprocessing of stress resultants)
  std::cerr <<" *** KirchhoffLoveShell::evalSol not implemented."<< std::endl;
  return false;
}


std::string KirchhoffLoveShell::getField1Name (size_t i,
                                               const char* prefix) const
{
  if (i >= 3) return "";

  char name = 'u'+i;
  if (!prefix)
    return std::string(1,name);

  return prefix + std::string(" ") + std::string(1,name);
}


std::string KirchhoffLoveShell::getField2Name (size_t i,
                                               const char* prefix) const
{
  if (i >= 6) return "";

  static const char* s[6] = { "n_xx", "n_yy", "n_xy", "m_xx", "m_yy", "m_xy" };

  if (!prefix)
    return s[i];

  return prefix + std::string(" ") + s[i];
}
