// $Id$
//==============================================================================
//!
//! \file SIMElasticity.C
//!
//! \date Dec 04 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#include "ElasticityUtils.h"
#include "SIMgeneric.h"
#include "IFEM.h"


//! Plane strain/stress option for 2D problems.
bool Elastic::planeStrain = false;
//! Axisymmtry option for 2D problems.
bool Elastic::axiSymmetry = false;
//! Option for Gauss point output to VTF.
bool Elastic::GIpointsVTF = false;
//! Time for function evaluation in linear problems.
double Elastic::time = 1.0;


void Elastic::printNorms (const Vector& gNorm, const Vector& rNorm,
                          const std::string& name, const SIMbase* model)
{
  const char* uRef = model->haveAnaSol() ? "|u|)   : " : "|u^r|) : ";
  double Rel = 100.0/(model->haveAnaSol() ? rNorm(3) : gNorm(1));

  IFEM::cout <<"\n\n>>> Error estimates based on "<< name <<" <<<"
             <<"\nEnergy norm |u^r| = a(u^r,u^r)^0.5   : "<< gNorm(1)
             <<"\nError norm a(e,e)^0.5, e=u^r-u^h     : "<< gNorm(2)
             <<"\n- relative error (% of "<< uRef << gNorm(2)*Rel;

  bool haveResErr = gNorm.size() >= (model->haveAnaSol() ? 8 : 5);
  if (haveResErr)
  {
    IFEM::cout <<"\nResidual error (r(u^r) + J(u^r))^0.5 : "<< gNorm(5)
               <<"\n- relative error (% of "<< uRef << gNorm(5)*Rel;
    if (gNorm.size() >= 9)
      IFEM::cout <<"\nJump term J(u^r)^0.5          : "<< gNorm(6)
                 <<"\n- relative error (% of "<< uRef << gNorm(6)*Rel;
  }

  if (model->haveAnaSol())
  {
    double exaErr = gNorm(gNorm.size() - (haveResErr ? 2 : 1));
    IFEM::cout <<"\nExact error a(e,e)^0.5, e=u-u^r      : "<< exaErr
               <<"\n- relative error (% of "<< uRef << exaErr*Rel
               <<"\nEffectivity index             : "<< gNorm(2)/rNorm(4);
    if (haveResErr)
      IFEM::cout <<"\nEffectivity index, theta^EX          : "
                 << (gNorm(2)+exaErr)/rNorm(4)
                 <<"\nEffectivity index, theta^RES         : "
                 << (gNorm(2)+gNorm(5))/rNorm(4);
  }
}


void Elastic::printBoundaryForces (const Vector& sf, RealArray& weights,
                                   const std::map<int,Vec3>& bCode,
                                   const SIMgeneric* model,
                                   bool indented)
{
  if (bCode.size() > 1 && weights.empty())
  {
    std::vector<int> glbNodes;
    weights.resize(model->getNoNodes());
    for (const std::pair<const int,Vec3>& c : bCode)
    {
      model->getBoundaryNodes(c.first,glbNodes);
      for (int inod : glbNodes) ++weights[inod-1];
    }
    for (double& w : weights)
      if (w > 1.0) w = 1.0/w;
  }

  int isec = 0;
  for (const std::pair<const int,Vec3>& c : bCode)
  {
    ++isec;
    Vector force = model->getInterfaceForces(sf,weights,c.first);
    if (force.normInf() > utl::zero_print_tol)
    {
      IFEM::cout << (indented ? "\n  " : "\n")
                 <<"Interface force at section "<< isec <<":";
      for (double f : force) IFEM::cout <<" "<< utl::trunc(f);
    }
  }

  IFEM::cout << std::endl;
}
