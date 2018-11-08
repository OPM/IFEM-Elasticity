// $Id$
//==============================================================================
//!
//! \file SIMElastic1D.C
//!
//! \date Nov 10 2018
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for 1D elastic solution drivers.
//!
//==============================================================================

#include "SIMElastic1D.h"
#include "SIMenums.h"
#include "ASMbase.h"
#include "AlgEqSystem.h"
#include "IntegrandBase.h"
#include "Vec3Oper.h"
#include "IFEM.h"


int SIMElastic1D::findLoadPoint (int ipt, int patch, double& u, bool onElement)
{
  IFEM::cout <<"\nLoad point #"<< ipt <<": patch #"<< patch <<" u="<< u;

  Vec3 X;
  int iclose = 0;
  int imatch = this->evalPoint(&u,X,&u,patch,true);
  if (imatch >= 0 && onElement)
    ipt = -this->findElementContaining(&u,patch);
  else if (imatch == 0 && (iclose = this->findClosestNode(X)) > 0)
    ipt = iclose;
  else
    ipt = imatch;

  if (ipt < 0 && onElement)
    IFEM::cout <<" on element #"<< -ipt <<" (u="<< u <<")";
  else if (iclose > 0)
    IFEM::cout <<", (closest) node #"<< ipt;
  else if (imatch > 0)
    IFEM::cout <<", node #"<< ipt;
  else
  {
    IFEM::cout << std::endl;
    std::cerr <<" *** Not a nodal point."<< std::endl;
    return 0;
  }

  IFEM::cout <<", X = "<< X;
  if (iclose > 0)
    IFEM::cout <<" (Xnod = "<< this->getNodeCoord(iclose) <<")";

  return ipt;
}


bool SIMElastic1D::assemblePoint (int patch, double u, double f, int ldof)
{
  ASMbase* pch = this->getPatch(patch,true);
  if (!pch) return false;

  Vec3 Fvec; Fvec(ldof) = f;
  this->setMode(SIM::RHS_ONLY);
  return pch->diracPoint(*myProblem,*myEqSys,&u,Fvec);
}
