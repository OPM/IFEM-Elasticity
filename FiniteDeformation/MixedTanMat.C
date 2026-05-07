// $Id$
//==============================================================================
//!
//! \file MixedTanMat.C
//!
//! \date Apr 29 2026
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Utilities for mixed model constitutive tangent matrix.
//!
//==============================================================================

#include "MixedTanMat.h"


/*!
  This method projects a 6x6 constitutive matrix into an equivalent
  7x7 constitutive matrix for the mixed model.
  \code
                | D_11  D_12 |
            D = |            |
                | D_21  D_22 |
  \endcode
  where:
  - D_11  = 6x6 Deviatoric part of matrix
  - D_12  = 6x1 Coupling   part of matrix
  - D_21  = 1x6 Coupling   part of matrix
  - D_22  = 1x1 Volumetric part of matrix
*/

void MixedMat::addVolumetricTerms (utl::matrix<double>& D,
                                   const std::vector<double>& sigma,
                                   double p_bar, double p_mix)
{
  size_t i, j;

  // Extend the matrix to 7x7 while keeping its content

  D.expandRows(7-D.rows());
  D.resize(7,7);

  // Compute left and right multiples with trace

  for (i = 1; i < 7; i++)
  {
    D(i,7) = (D(i,1) + D(i,2) + D(i,3)) / 3.0;
    D(7,i) = (D(1,i) + D(2,i) + D(3,i)) / 3.0;
  }

  // Convert upper 6x6 matrix to a the deviatoric D_11

  for (i = 1; i < 7; i++)
    for (j = 1; j <= 3; j++)
    {
      D(i,j) -= D(i,7);
      D(j,i) -= D(7,i);
    }

  // Form the last term, D_22

  D(7,7) = (D(1,7) + D(2,7) + D(3,7)) / 3.0;

  // Final update to form D_12 and D_21

  for (i = 1; i <= 3; i++)
  {
    D(i,7) -= D(7,7);
    D(7,i) -= D(7,7);
    for (j = 1; j <= 3; j++)
      D(i,j) += D(7,7);
  }

  // Compute deviatoric stress

  utl::vector<double> sig_dev(sigma);
  for (size_t i = 1; i <= 3; i++)
    sig_dev(i) -= p_bar;
  for (double& sd : sig_dev)
    sd *= 2.0/3.0;

  // D_11: B_matrix part

  const double f1 = p_mix - p_bar*2.0/3.0;
  const double f2 = p_bar - p_mix;

  for (i = 1; i <= 3; i++)
    for (j = 1; j <= 3; j++)
      D(i,j) += f1;

  for (i = 1; i <= 3; i++)
    for (j = 1; j <= 6; j++)
    {
      D(i,j) -= sig_dev(j);
      D(j,i) -= sig_dev(j);
    }

  for (j = 1; j <= 3; j++)
  {
    D(j  ,j  ) += 2.0*f2;
    D(j+3,j+3) += f2;
  }

  // D_12: Coupling matrix with

  for (j = 1; j < 7; j++)
  {
    D(7,j) += sig_dev(j);
    D(j,7) += sig_dev(j);
  }

  // D_22: Volumetric part

  D(7,7) -= p_bar/3.0;
}
