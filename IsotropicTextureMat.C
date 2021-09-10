// $Id$
//==============================================================================
//!
//! \file IsotropicTextureMat.C
//!
//! \date Sep 13 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Isotropic linear elastic material model defined through a texture.
//!
//==============================================================================

#include "IsotropicTextureMat.h"
#include "FiniteElement.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"
#include "StbImage.h"


void IsotropicTextureMat::parse (const TiXmlElement* elem)
{
  std::string textureFile;
  utl::getAttribute(elem, "file", textureFile);

  int width, height, nrChannels;
  unsigned char* image = stb::loadImage(textureFile.c_str(),
                                        width, height, nrChannels);
  if (!image) {
    std::cerr << "File not found: " << textureFile << std::endl;
    return;
  }

  textureData.resize(width,std::vector<rgba>(height,{0.0,0.0,0.0}));
  const unsigned char* data = image;
  for (int j = 0; j < height; j++)
    for (int i = 0; i < width; i++)
      for (int c = 0; c < nrChannels; c++)
        textureData[i][j][c] = double(*data++) / 255.0;

  free(image);

  Doubles     range;
  LinIsotropic mat(planeStress,axiSymmetry);
  const TiXmlElement* child = elem->FirstChildElement("range");
  for (; child; child = child->NextSiblingElement("range"))
  {
    utl::getAttribute(child,"min",range.first);
    utl::getAttribute(child,"max",range.second);
    IFEM::cout << (materials.empty() ? "\n\t" : "\t");
    mat.parse(child);
    materials[range] = mat;
  }
}


void IsotropicTextureMat::printLog () const
{
  for (const std::pair<const Doubles,LinIsotropic>& mat : materials)
  {
    IFEM::cout <<"Material with range ["
               << mat.first.first <<","<< mat.first.second <<"]:\n";
    mat.second.printLog();
  }
}


const LinIsotropic* IsotropicTextureMat::findMaterial (const FiniteElement& fe) const
{
  if (textureData.empty())
    return nullptr;

  int nrow = textureData.size();
  int ncol = textureData.front().size();
  int i = fe.u * (nrow-1);
  int j = fe.v * (ncol-1);
  if (i < 0 || i >= nrow || j < 0 || j >= ncol)
  {
    std::cerr <<" *** Texture index out of bounds "<< i <<" "<< j << std::endl;
    return nullptr;
  }

  double I = textureData[i][j].front();
  auto mat = std::find_if(materials.begin(),materials.end(),
                          [I](const std::pair<Doubles,LinIsotropic>& a)
                          {
                            return a.first.first <= I && I <= a.first.second;
                          });

  return mat == materials.end() ? nullptr : &mat->second;
}


bool IsotropicTextureMat::evaluate (double& lambda, double& mu,
                                    const FiniteElement& fe,
                                    const Vec3& X) const
{
  const LinIsotropic* mat = this->findMaterial(fe);
  return mat ? mat->evaluate(lambda, mu, fe, X) : false;
}


bool IsotropicTextureMat::evaluate (Matrix& C, SymmTensor& sigma, double& U,
                                    const FiniteElement& fe, const Vec3& X,
                                    const Tensor&, const SymmTensor& eps,
                                    char iop, const TimeDomain*, const Tensor*) const
{
  const LinIsotropic* mat = this->findMaterial(fe);
  return mat ? mat->evaluate(C, sigma, U, fe, X, eps, eps, iop) : false;
}
