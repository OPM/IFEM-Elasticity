// $Id$
//==============================================================================
//!
//! \file SIMElasticityWrap.C
//!
//! \date April 21 2016
//!
//! \author Knut Morten Okstad
//!
//! \brief Wrapper equipping the elasticity solver with an ISolver interface.
//!
//==============================================================================

#include "SIMElasticityWrap.h"

#include "DataExporter.h"
#include "SIM2D.h"
#include "SIM3D.h"


template<class Dim>
SIMElasticityWrap<Dim>::SIMElasticityWrap ()
{
  Dim::msgLevel = 1;
  Dim::myHeading = "Elasticity solver";
}


template<class Dim>
void SIMElasticityWrap<Dim>::registerFields (DataExporter& exporter)
{
  int flag = DataExporter::PRIMARY;
  if (!Dim::opt.pSolOnly)
    flag |= DataExporter::SECONDARY;
  exporter.registerField("u","solution",DataExporter::SIM,flag);
  exporter.setFieldValue("u",this,&this->getSolution());
}


template<class Dim>
bool SIMElasticityWrap<Dim>::saveModel (char* fileName, int& geoBlk, int&)
{
  if (Dim::opt.format < 0)
    return true;

  return this->writeGlvG(geoBlk,fileName);
}


template<class Dim>
bool SIMElasticityWrap<Dim>::saveStep (const TimeStep& tp, int& nBlock)
{
  if (Dim::opt.format < 0 || tp.step%Dim::opt.saveInc > 0)
    return true;

  int iDump = 1 + tp.step/Dim::opt.saveInc;
  if (!this->writeGlvS(this->getSolution(),iDump,nBlock,tp.time.t,"u"))
    return false;

  return this->writeGlvStep(iDump,tp.time.t);
}


template<class Dim>
bool SIMElasticityWrap<Dim>::serialize (SerializeMap& data) const
{
  if (!this->saveBasis(data) || !this->saveSolution(data,this->getName()))
    return false;

  data["Elasticity::Eext"] = SIMsolution::serialize(this->getExtEnerg(),1);

  return true;
}


template<class Dim>
bool SIMElasticityWrap<Dim>::deSerialize (const SerializeMap& data)
{
  if (!this->restoreSolution(data,this->getName()))
    return false;

  SerializeMap::const_iterator sit = data.find("Elasticity::Eext");
  if (sit != data.end())
    SIMsolution::deSerialize(sit->second,this->theExtEnerg(),1);

  return true;
}


template<class Dim>
bool SIMElasticityWrap<Dim>::deSerializeBasis (const SerializeMap& data)
{
  return this->restoreBasis(data);
}


template<class Dim>
bool SIMElasticityWrap<Dim>::init (const TimeStep& tp, bool withRF)
{
  return (this->initSystem(Dim::opt.solver,1,1,0,withRF) &&
          this->initSolution(this->getNoDOFs(),this->getNoSolutions()) &&
          this->setMode(SIM::INIT) &&
          this->getIntegrand()->init(tp.time));
}


template<class Dim>
bool SIMElasticityWrap<Dim>::advanceStep (TimeStep& tp)
{
  this->pushSolution(); // Update solution vectors between time steps
  return this->SIMElasticity<Dim>::advanceStep(tp);
}


template class SIMElasticityWrap<SIM2D>;
template class SIMElasticityWrap<SIM3D>;
