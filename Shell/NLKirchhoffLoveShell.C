// $Id$
//==============================================================================
//!
//! \file NLKirchhoffLoveShell.C
//!
//! \date Apr 20 2018
//!
//! \author Simen Skogholt Haave and Marit Gaarder Rakvaag / NTNU
//!
//! \brief Class for nonlinear Kirchhoff-Love thin shell problems.
//!
//==============================================================================

#include "NLKirchhoffLoveShell.h"
#include "FiniteElement.h"
#include "ElmMats.h"
#include "CoordinateMapping.h"
#include "Vec3Oper.h"


void NLKirchhoffLoveShell::setMode (SIM::SolutionMode mode)
{
  this->KirchhoffLoveShell::setMode(mode);
  if (mode == SIM::STATIC || mode == SIM::RHS_ONLY)
  {
    iS = 1;
    primsol.resize(1);
  }
}


int NLKirchhoffLoveShell::getIntegrandType () const
{
  return SECOND_DERIVATIVES | UPDATED_NODES;
}


bool NLKirchhoffLoveShell::evalInt (LocalIntegral& elmInt,
                                    const FiniteElement& fe,
                                    const Vec3& X) const
{
  Matrix Gd, Hd;
  if (elmInt.vec.size() > 1)
  {
    // Co-variant basis and Hessian in deformed configuration
    Matrix3D Hess;
    Gd.multiplyMat(elmInt.vec.back(),fe.dNdX);
    if (Hess.multiplyMat(elmInt.vec.back(),fe.d2NdX2))
      utl::Hessian(Hess,Hd);
    else
      return false;
  }

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (eM) // Integrate the mass matrix
    this->formMassMatrix(elMat.A[eM-1],fe.N,X,fe.detJxW);

  if (eK && iS) // Integrate the stiffness matrix and internal forces
    this->evalKandS(elMat.A[eK-1],elMat.b[iS-1],fe,fe.G,Gd,fe.H,Hd,X);
  else if (iS) // Integrate the internal forces only
  {
    Matrix dummyEK;
    this->evalKandS(dummyEK,elMat.b[iS-1],fe,fe.G,Gd,fe.H,Hd,X);
  }

  if (eS) // Integrate the load vector due to gravitation and other body forces
  {
    Vec3 n;
    if (presFld)
      n = this->getShellNormal(Gd.empty() ? fe.G : Gd);
    this->formBodyForce(elMat.b[eS-1],fe.N,fe.iGP,X,n,fe.detJxW);
  }

  return true;
}


bool NLKirchhoffLoveShell::evalKandS (Matrix& EK, Vector& ES,
                                      const FiniteElement& fe,
                                      const Matrix& G0, const Matrix& Gn,
                                      const Matrix& H0, const Matrix& Hn,
                                      const Vec3& X) const
{
  Matrix Dm, Db;
  if (!this->formDmatrix(Dm,Db,fe,X))
    return false;

  return true;
}
