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
#include "IFEM.h"
#include "tinyxml.h"
#include "FiniteElement.h"
#include "Utilities.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

void IsotropicTextureMat::parse (const TiXmlElement* elem)
{
  std::string textureFile;
  utl::getAttribute(elem, "file", textureFile);

  int width, height, nrChannels;
  unsigned char *image = stbi_load(textureFile.c_str(),  &width, &height, &nrChannels, 0);
  if (!image) {
    std::cerr << "File not found: " << textureFile << std::endl;
    return;
  }
  textureData.resize(width, std::vector<std::array<double,4> >(height));
  const unsigned char* data = image;
  for(int j=0; j<height; ++j) {
    for(int i=0; i<width; ++i) {
      for(int c=0; c<nrChannels; c++) {
        textureData[i][j][c] = double(*data++) / 255.0;
      }
      for(int c=nrChannels; c<4; c++) {
        textureData[i][j][c] = textureData[i][j][c];
      }
    }
  }
  free(image);

  const TiXmlElement* child = elem->FirstChildElement("range");
  while (child) {
    std::pair<double,double> minmax;
    utl::getAttribute(child,"min",minmax.first);
    utl::getAttribute(child,"max",minmax.second);

    materials.insert(std::make_pair(minmax,LinIsotropic(!planeStress,axiSymmetry)));
    materials[minmax].parse(child);
    child = child->NextSiblingElement("range");
  }
  for (const auto it : materials) {
    IFEM::cout << "Material with range [" << it.first.first << "," << it.first.second << "]:\n";
    it.second.printLog();
  }
}


void IsotropicTextureMat::printLog () const
{
  for (const auto it : materials) {
    IFEM::cout << "Material with range [" << it.first.first << "," << it.first.second << "]:\n";
    it.second.printLog();
  }
}


const LinIsotropic* IsotropicTextureMat::findMaterial (const FiniteElement& fe) const
{
  double I = this->findIntensity(fe);

  auto mat  = std::find_if(materials.begin(), materials.end(),
            [I](const std::pair<std::pair<double,double>,LinIsotropic>& a)
            {
              return a.first.first <= I && I <= a.first.second;
            });

  if (mat == materials.end())
    return nullptr;

  return &mat->second;

}


double IsotropicTextureMat::findIntensity (const FiniteElement& fe) const
{
  size_t i = fe.u * (textureData.size()   -1);
  size_t j = fe.v * (textureData[0].size()-1);
  if (i<0 || i>textureData.size()   -1 ||
      j<0 || j>textureData[0].size()-1) {
    std::cerr << "Texture index out of bounds" << std::endl;
    return -1;
  }
  return textureData[i][j][0];
}


bool IsotropicTextureMat::evaluate (double& lambda, double& mu,
                                    const FiniteElement& fe,
                                    const Vec3& X) const
{
  const LinIsotropic* mat = this->findMaterial(fe);
  if (!mat)
    return false;

  return mat->evaluate(lambda, mu, fe, X);
}


bool IsotropicTextureMat::evaluate (Matrix& C, SymmTensor& sigma, double& U,
                                    const FiniteElement& fe, const Vec3& X,
                                    const Tensor& T, const SymmTensor& eps,
                                    char iop, const TimeDomain* time, const Tensor* T2) const
{
  const LinIsotropic* mat = this->findMaterial(fe);
  if (!mat)
    return false;

  return mat->evaluate(C, sigma, U, fe, X, T, eps, iop, time, T2);
}
