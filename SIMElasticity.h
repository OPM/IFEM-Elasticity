// $Id$
//==============================================================================
//!
//! \file SIMElasticity.h
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#ifndef _SIM_ELASTICITY_H
#define _SIM_ELASTICITY_H

#include "SIMRigid.h"
#include "MatVec.h"
#include "Vec3.h"

class ElasticBase;
class Material;
class TimeStep;
struct TimeDomain;

typedef std::vector<Material*> MaterialVec; //!< Convenience declaration


/*!
  \brief Driver class for isogeometric FEM analysis of elasticity problems.
  \details The class incapsulates data and methods for solving elasticity
  problems using NURBS-based finite elements. It reimplements the parse methods
  and some property initialization methods of the parent class.
*/

template<class Dim> class SIMElasticity : public Dim, protected SIMRigid
{
public:
  //! \brief Default constructor.
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  explicit SIMElasticity(bool checkRHS = false);

  //! \brief The destructor frees the dynamically allocated material properties.
  virtual ~SIMElasticity();

  //! \brief Returns the name of this simulator (for use in the HDF5 export).
  virtual std::string getName() const { return "Elasticity"; }

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp);

  //! \brief Initializes the property containers of the model.
  virtual void clearProperties();

  //! \brief Prints out load step identification.
  //! \param[in] istep Load step counter
  //! \param[in] time Parameters for nonlinear simulations
  virtual void printStep(int istep, const TimeDomain& time) const;

  //! \brief Calculates surface traction resultants.
  //! \param[out] f Calculated traction resultants
  //! \param[in] sol Primary solution vectors
  virtual bool calcBouForces(Vectors& f, const Vectors& sol);

  //! \brief Calculates the traction resultant associated with a given boundary.
  //! \param[out] f Calculated traction resultant
  //! \param[in] sol Primary solution vectors
  //! \param[in] tp Time stepping parameters
  bool getBoundaryForce(Vector& f, const Vectors& sol, const TimeStep& tp);
  //! \brief Extracts reaction forces associated with given boundaries.
  //! \param[out] rf Reaction force resultant for specified boundaries
  bool getBoundaryReactions(Vectors& rf);
  //! \brief Extracts reaction forces associated with given boundary.
  //! \param[out] rf Reaction force resultant for the specified boundary
  //! \param[in] bindex One-based boundary code index, zero for the sum
  bool getBoundaryReactions(Vector& rf, size_t bindex = 0);

  //! \brief Returns whether reaction forces are to be computed or not.
  virtual bool haveBoundaryReactions(bool reactionsOnly = false) const;
  //! \brief Returns whether an analytical solution is available or not.
  virtual bool haveAnaSol() const;

protected:
  //! \brief Performs some preprocessing tasks before the FEM model generation.
  virtual void preprocessA();
  //! \brief Specialized preprocessing performed before assembly initialization.
  virtual bool preprocessBeforeAsmInit(int& ngnod);
  //! \brief Preprocessing performed after the system assembly initialization.
  virtual bool preprocessB();

  //! \brief Returns the actual integrand.
  virtual ElasticBase* getIntegrand() = 0;

  //! \brief Parses the analytical solution from an input stream.
  virtual bool parseAnaSol(char*, std::istream&);
  //! \brief Parses the analytical solution from an XML element.
  virtual bool parseAnaSol(const tinyxml2::XMLElement*);

  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML element
  //! \param[in] elem The XML element to parse
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Preprocesses a user-defined Dirichlet boundary property.
  //! \param[in] patch 1-based index of the patch to receive the property
  //! \param[in] lndx Local index of the boundary item to receive the property
  //! \param[in] ldim Dimension of the boundary item to receive the property
  //! \param[in] dirs Which local DOFs to constrain
  //! \param[in] code In-homegeneous Dirichlet condition property code
  //! \param ngnod Total number of global nodes in the model (might be updated)
  //! \param[in] basis Which basis to apply the constraint to (mixed methods)
  virtual bool addConstraint(int patch, int lndx, int ldim,
                             int dirs, int code, int& ngnod, char basis);

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd);

  //! \brief Initializes the body load properties for current patch.
  //! \param[in] patchInd 1-based patch index
  virtual bool initBodyLoad(size_t patchInd);

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd);

  //! \brief Returns norm index of the integrated volume.
  virtual size_t getVolumeIndex() const;

  //! \brief Reverts the square-root operation on the volume and VCP quantities.
  virtual bool postProcessNorms(Vectors& gNorm, Matrix* eNorm);

public:
  //! \brief Prints a norm group to the log stream.
  //! \param[in] gNorm The norm values to print
  //! \param[in] rNorm Reference norms for the first norm group
  //! \param[in] prjName Projection name associated with this norm group
  virtual void printNormGroup(const Vector& gNorm, const Vector& rNorm,
                              const std::string& prjName) const;

  //! \brief Prints interface force resultants associated with given boundaries.
  //! \param[in] sf Internal nodal forces
  //! \param weights Nodal weights (in case some nodes are present in more sets)
  virtual void printIFforces(const Vector& sf, RealArray& weights);

  //! \brief Writes current model geometry to the VTF-file.
  //! \param nBlock Running result block counter
  //! \param[in] inpFile File name used to construct the VTF-file name from
  //! \param[in] doClear If \e true, clear geometry block if \a inpFile is null
  virtual bool writeGlvG(int& nBlock, const char* inpFile, bool doClear = true);

protected:
  MaterialVec mVec;      //!< Material data
  std::string myContext; //!< XML-tag to search for problem inputs within
  std::map<int,Vec3> bCode; //!< Property codes for boundary traction resultants

private:
  bool plotRgd; //!< If \e true, output rigid couplings as VTF geometry
  int  aCode;   //!< Analytical BC code (used by destructor)
};

#endif
