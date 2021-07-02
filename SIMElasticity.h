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


class Elasticity;
class Material;
class TimeStep;


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

  //! \brief Calculates surface traction resultants.
  //! \param[out] f Calculated traction resultants
  //! \param[in] sol Primary solution vectors
  //!
  //! \details The boundaries for which the traction is calculated are
  //! identified by the property set codes in \a bCode, which are
  //! assigned values by parsing `<boundaryforce>` tags in the input file.
  virtual bool calcBouForces(Vectors& f, const Vectors& sol);

  //! \brief Calculates the traction resultant associated with a given boundary.
  //! \param[out] f Calculated traction resultant
  //! \param[in] sol Primary solution vectors
  //! \param[in] tp Time stepping parameters
  //!
  //! \details The boundary for which the traction is calculated is identified
  //! by the property set code \a bCode which is assigned value by parsing
  //! the first `<boundaryforce>` tag in the input file.
  bool getBoundaryForce(Vector& f, const Vectors& sol, const TimeStep& tp);

  //! \brief Extracts the reaction forces associated with a given boundary.
  //! \param[out] rf Reaction force resultant for specified boundary
  //!
  //! \details The boundary for which the reaction force is returned
  //! is identified by the property set code \a bCode which is assigned value
  //! by parsing the first `<boundaryforce>` tag in the input file.
  bool getBoundaryReactions(Vector& rf);

  //! \brief Returns whether reaction forces are to be computed or not.
  bool haveBoundaryReactions() const;

  //! \brief Returns whether an analytical solution is available or not.
  virtual bool haveAnaSol() const;

protected:
  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented inserting a call to getIntegrand().
  //! This makes sure the integrand has been allocated in case of minimum input.
  //! It also resolves inhomogeneous boundary condition fields in case they are
  //! derived from the analytical solution.
  virtual void preprocessA();

  //! \brief Specialized preprocessing performed before assembly initialization.
  //! \details This method creates the multi-point constraint equations
  //! representing the rigid couplings in the model.
  virtual bool preprocessBeforeAsmInit(int& ngnod);

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to ensure that threading groups are
  //! established for the patch faces subjected to boundary force integration.
  //! In addition, the reference point for moment calculation \b X0 of each
  //! boundary is calculated based on the control/nodal point coordinates.
  virtual bool preprocessB();

  //! \brief Returns the actual integrand.
  virtual Elasticity* getIntegrand() = 0;

  //! \brief Parses the analytical solution from an input stream.
  virtual bool parseAnaSol(char*, std::istream&);

  //! \brief Parses the analytical solution from an XML element.
  virtual bool parseAnaSol(const TiXmlElement*);

  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Parses a data section from an XML element
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Preprocesses a user-defined Dirichlet boundary property.
  //! \param[in] patch 1-based index of the patch to receive the property
  //! \param[in] lndx Local index of the boundary item to receive the property
  //! \param[in] ldim Dimension of the boundary item to receive the property
  //! \param[in] dirs Which local DOFs to constrain
  //! \param[in] code In-homegeneous Dirichlet condition property code
  //! \param ngnod Total number of global nodes in the model (might be updated)
  //! \param[in] basis Which basis to apply the constraint to (mixed methods)
  //!
  //! \details This method is overridden to handle dirichlet conditions on
  //! the explicit master nodes of rigid couplings which not are regular nodes
  //! in a patch. These nodes may also have rotational degrees of freedom.
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
  virtual void printNormGroup (const Vector& gNorm, const Vector& rNorm,
                               const std::string& prjName) const;

  //! \brief Writes current model geometry to the VTF-file.
  //! \param nBlock Running result block counter
  //! \param[in] inpFile File name used to construct the VTF-file name from
  //! \param[in] doClear If \e true, clear geometry block if \a inpFile is null
  //!
  //! \details This method is overrriden to also account for rigid couplings.
  virtual bool writeGlvG(int& nBlock, const char* inpFile, bool doClear = true);

protected:
  MaterialVec mVec;      //!< Material data
  std::string myContext; //!< XML-tag to search for problem inputs within

private:
  bool plotRgd; //!< If \e true, output rigid couplings as VTF geometry
  int  aCode;   //!< Analytical BC code (used by destructor)
  std::map<int,Vec3> bCode; //!< Property codes for boundary traction resultants
};

#endif
