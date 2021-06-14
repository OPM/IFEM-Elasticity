// $Id$
//==============================================================================
//!
//! \file TestStaticCondensation.C
//!
//! \date May 6 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Unit tests for static condensation of linear systems.
//!
//==============================================================================

#include "SIMLinEl.h"
#include "DenseMatrix.h"
#include "SAM.h"
#include <array>

#include "gtest/gtest.h"


TEST(TestSIMLinEl2D, StaticCondensation)
{
  SIMLinEl2D fullModel(false,false), scModel(false,false);
  Vectors    displ(2);

  // Initialize the full model from input file
  ASSERT_TRUE(fullModel.read("SSmembrane-p1.xinp"));
  ASSERT_TRUE(fullModel.preprocess());

  // Assemble and direct solve the full equation system
  fullModel.setMode(SIM::STATIC);
  fullModel.setQuadratureRule(fullModel.opt.nGauss[0],true,true);
  fullModel.initSystem(LinAlg::SPARSE);
  ASSERT_TRUE(fullModel.assembleSystem());
  ASSERT_TRUE(fullModel.solveSystem(displ,1));
  std::cout <<"Full Solution:"<< displ.front();

  // Initialize the model to be solved by static condensation
  SIMbase::ignoreDirichlet = true;
  ASSERT_TRUE(scModel.read("SSmembrane-p1.xinp"));
  ASSERT_TRUE(scModel.preprocess());

  // Assemble equation system and perform static condensation
  Matrix Kred;
  Vector Rred;
  scModel.opt.num_threads_SLU = -1; // To avoid pre-assembly
  ASSERT_TRUE(scModel.staticCondensation(Kred,Rred));
#ifdef INT_DEBUG
  std::cout <<"Kred:"<< Kred;
  std::cout <<"Rred:"<< Rred;
#endif

  // Enforce boundary conditions (super DOFs 1, 2 and 5)
  std::array<size_t,3> fixed{1,2,5};
  for (int i : fixed)
  {
    for (size_t j = 1; j <= Kred.cols(); j++)
      Kred(i,j) = Kred(j,i) = 0.0;
    Kred(i,i) = 1.0;
    Rred(i) = 0.0;
  }
#ifdef INT_DEBUG
  std::cout <<"Constrained Kred:"<< Kred;
  std::cout <<"Constrained Rred:"<< Rred;
#endif

  // Solve the condensed system
  DenseMatrix Ksup(Kred);
  StdVector   Rsup(Rred);
  ASSERT_TRUE(Ksup.solve(Rsup,1));
  std::cout <<"Superelement Solution:"<< Rsup;

  // Find global DOF numbers of the "supernodes" 106 and 107
  std::pair<int,int> dofs1 = fullModel.getSAM()->getNodeDOFs(106);
  std::pair<int,int> dofs2 = fullModel.getSAM()->getNodeDOFs(107);
  std::cout <<"Node 106: "<< dofs1.first <<" - "<< dofs1.second << std::endl;
  std::cout <<"Node 107: "<< dofs2.first <<" - "<< dofs2.second << std::endl;

  // Compare the solutions, they should be identical
  int idof, i = 0;
  for (idof = dofs1.first; idof <= dofs1.second; idof++, i++)
    EXPECT_NEAR(displ.front()(idof), Rsup[i], 1.0e-8);
  for (idof = dofs2.first; idof <= dofs2.second; idof++, i++)
    EXPECT_NEAR(displ.front()(idof), Rsup[i], 1.0e-8);

  // Check that the recovery matrix is valid
  FILE* fd = fopen("SSmembrane-p1.rec","rb");
  ASSERT_FALSE(fd == nullptr);

  size_t len = 0, n1, n2;
  char* header = nullptr;
  ASSERT_GT(getline(&header,&len,fd),0);
  ASSERT_FALSE(strstr(header,"#IFEM recovery matrix:") == nullptr);
  EXPECT_EQ(sscanf(header+22,"%zu%zu",&n1,&n2),2);
  std::cout << header << n1 <<" "<< n2;
  free(header);

  Matrix Rmat(n1,n2);
  for (size_t c = 0; c < n2; c++)
    EXPECT_EQ(fread(Rmat.ptr(c),sizeof(double),n1,fd),n1);
  fclose(fd);
#if INT_DEBUG > 2
  std::cout << Rmat;
#else
  std::cout << std::endl;
#endif

  // Recover internal displacements and verify against the full solution
  ASSERT_TRUE(scModel.recoverInternals(Rsup,displ.back()));
  std::cout <<"Recovered Solution:"<< displ.back();
  ASSERT_EQ(displ.front().size(),displ.back().size());
  for (size_t i = 0; i < displ.front().size(); i++)
    EXPECT_NEAR(displ.front()[i], displ.back()[i], 1.0e-8);
}
