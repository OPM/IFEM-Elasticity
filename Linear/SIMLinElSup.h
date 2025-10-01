// $Id$
//==============================================================================
//!
//! \file SIMLinElSup.h
//!
//! \date Jun 12 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for linear elastic superelement FEM analysis.
//!
//==============================================================================

#ifndef _SIM_LIN_EL_SUP_H
#define _SIM_LIN_EL_SUP_H

#include "SIMsupel.h"


/*!
  \brief Solution driver for linear elastic superelement FEM analysis.
*/

class SIMLinElSup : public SIMsupel
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  //! \param[in] hd Sub-simulator heading
  //! \param[in] fd If \e true, merge duplicated FE nodes on patch interfaces
  SIMLinElSup(const char* hd, bool fd) : SIMsupel(hd), fixDup(fd) {}
  //! \brief The destructor deletes the FE substructure data.
  virtual ~SIMLinElSup();

protected:
  using SIMsupel::parse;
  //! \brief Parses a data section from an XML element
  //! \param[in] elem The XML element to parse
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to preprocess the FE substructures.
  virtual bool preprocessB();

  //! \brief Tesselates the superelement associated with specified patch.
  virtual ElementBlock* tesselatePatch(size_t pidx) const;

  //! \brief Returns a pointer to the problem-specific data object.
  IntegrandBase* getMyProblem() const;

  //! \brief Writes primary solution for a given load case to the VTF-file.
  //! \param[in] psol Primary solution vector
  //! \param[in] iStep Load case identifier
  //! \param nBlock Running result block counter
  //! \param[in] idBlock Starting value of result block numbering
  virtual int writeGlvS1(const Vector& psol, int iStep, int& nBlock,
                         double, const char*, int idBlock, int, bool);

  //! \brief Writes secondary solution for a given load case to the VTF-file.
  //! \param[in] iStep Load case identifier
  //! \param nBlock Running result block counter
  //! \param[in] idBlock Starting value of result block numbering
  virtual int writeGlvS2(const Vector&, int iStep, int& nBlock,
                         double, int idBlock, int);

private:
  //! \brief Struct representing a FE substructure.
  struct FEmodel
  {
    SIMgeneric*   sim = nullptr; //!< Underlying FE model of the substructure
    ElementBlock* blk = nullptr; //!< Tesselated FE model for visualization
  };

  std::map<std::string,FEmodel> mySubSim; //!< FE substructure container

  bool fixDup; //!< If \e true, merge duplicated FE nodes on patch interfaces
};

#endif
