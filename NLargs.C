// $Id$
//==============================================================================
//!
//! \file NLargs.C
//!
//! \date Nov 20 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Pre-parsing of input files for Finite Deformation applications.
//!
//==============================================================================

#include "NLargs.h"
#include "NLoptions.h"
#include "SIMenums.h"
#include "ASMmxBase.h"
#include "ElasticityUtils.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml2.h"


bool NLargs::parseArg (const char* argv)
{
  int form = -1;
  int pOrd = -1;

  if (!strcmp(argv,"-printMax"))
    printMax = 'G';
  else if (!strcmp(argv,"-printMaxPatch"))
    printMax = 'P';
  else if (!strncmp(argv,"-dumpNod",8))
    dNodeMap = true;
  else if (!strcmp(argv,"-checkRHS"))
    checkRHS = true;
  else if (!strcmp(argv,"-fixDup"))
    fixDup = true;
  else if (!strncmp(argv,"-2Dpstra",8))
    twoD = Elastic::planeStrain = true;
  else if (!strncmp(argv,"-2Daxi",6))
    twoD = Elastic::axiSymmetry = true;
  else if (!strncmp(argv,"-2D",3))
    twoD = true;
  else if (!strcmp(argv,"-linear"))
    form = SIM::LINEAR;
  else if (!strcmp(argv,"-UL"))
  {
    if (form < SIM::UPDATED_LAGRANGE)
      form = SIM::UPDATED_LAGRANGE;
  }
  else if (!strncmp(argv,"-MX",3))
  {
    form = SIM::MIXED_QnPn1;
    if (strlen(argv) > 3 && isdigit(argv[3]))
      pOrd = atoi(argv+3);
  }
  else if (!strcmp(argv,"-Mixed"))
  {
    form = SIM::MIXED_QnQn1;
    ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;
  }
  else if (!strcmp(argv,"-mixed"))
  {
    form = SIM::MIXED_QnQn1;
    ASMmxBase::Type = ASMmxBase::REDUCED_CONT_RAISE_BASIS1;
  }
  else if (!strncmp(argv,"-Fbar",5))
  {
    form = SIM::FBAR;
    if (strlen(argv) > 5 && isdigit(argv[5]))
      pOrd = atoi(argv+5);
  }
  else if (!strcmp(argv,"-Tensor"))
    form = SIM::NONLINEAR;
  else if (!strcmp(argv,"-GA"))
    algor = GENALPHA;
  else if (!strcmp(argv,"-oldHHT"))
    algor = OLDHHT;
  else if (!strcmp(argv,"-HHT"))
    algor = NEWHHT;
  else if (!strcmp(argv,"-arclen"))
    algor = ARCLEN;
  else if (!strncmp(argv,"-adap",5))
    adaptive = true;
  else
    return false;

  this->setFormulation(form,pOrd);
  return true;
}


bool NLargs::parse (const tinyxml2::XMLElement* elem)
{
  if (!strcasecmp(elem->Value(),"geometry"))
  {
    int dim = 3;
    if (!utl::getAttribute(elem,"dimension",dim))
      utl::getAttribute(elem,"dim",dim);
    twoD = dim == 2;
  }
  else if (!strcasecmp(elem->Value(),"finitedeformation"))
  {
    utl::getAttribute(elem,"adaptive",adaptive);
    const tinyxml2::XMLElement* child = elem->FirstChildElement("formulation");
    if (child) this->parseFormulation(child);
  }
  else if (!strcasecmp(elem->Value(),"nonlinearsolver"))
  {
    if (elem->FirstChildElement("arclen"))
      algor = ARCLEN;
  }
  else if (!strcasecmp(elem->Value(),"newmarksolver"))
  {
    std::string version("hht");
    utl::getAttribute(elem,"version",version);
    if (version.find("alpha") != std::string::npos)
      algor = GENALPHA;
    else if (version.find("old") != std::string::npos)
      algor = OLDHHT;
    else
      algor = NEWHHT;
  }

  return IFEM::getOptions().parseDiscretizationTag(elem);
}


void NLargs::parseFormulation (const tinyxml2::XMLElement* elem)
{
  int form = -1;
  int pOrd = -1;

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(),"planestrain"))
      Elastic::planeStrain = twoD;
    else if (!strcasecmp(child->Value(),"axisymmetric"))
      Elastic::axiSymmetry = twoD;
    else if (!strcasecmp(child->Value(),"totallagrange"))
      form = SIM::TOTAL_LAGRANGE;
    else if (!strcasecmp(child->Value(),"updatedlagrange") &&
             form < SIM::UPDATED_LAGRANGE)
      form = SIM::UPDATED_LAGRANGE;
    else if (!strcasecmp(child->Value(),"linear"))
      form = SIM::LINEAR;
    else if (!strcasecmp(child->Value(),"mixed"))
    {
      if (child->FirstChild())
        pOrd = atoi(child->FirstChild()->Value());
      std::string type;
      utl::getAttribute(child,"type",type);
      if (type == "Qp/Pp-1")
        form = SIM::MIXED_QnPn1;
      else if (type == "Qp/Qp-1")
      {
        form = SIM::MIXED_QnQn1;
        if (pOrd == 1)
          ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;
        else
          ASMmxBase::Type = ASMmxBase::REDUCED_CONT_RAISE_BASIS1;
      }
      else if (type == "Fbar")
        form = SIM::FBAR;
      else if (type == "Tensor")
        form = SIM::NONLINEAR;
    }

  this->setFormulation(form,pOrd);
}


void NLargs::setFormulation (int form, int pOrd)
{
  if (form >= 0)
  {
    if (options.empty())
      options = { form };
    else
      options.front() = form;
  }

  if (pOrd >= 0)
  {
    if (options.empty())
      options = { 0, pOrd };
    else if (options.size() < 2)
      options.resize(2,pOrd);
    else
      options[1] = pOrd;
  }
}
