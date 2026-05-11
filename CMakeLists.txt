cmake_minimum_required(VERSION 3.22)

project(LinEl LANGUAGES C CXX)

message(STATUS "Configuring project ${PROJECT_NAME}")

if(NOT TARGET IFEM)
  if(IFEM_AS_SUBMODULE)
    add_subdirectory(../.. IFEM)
    set(IFEM_PATH ${PROJECT_SOURCE_DIR}/../..)
  else()
    get_filename_component(BUILD_DIR ${PROJECT_BINARY_DIR} NAME)
    set(IFEM_DIR ${PROJECT_SOURCE_DIR}/../../${BUILD_DIR})
    find_package(IFEM REQUIRED)
  endif()
endif()

if(NOT TARGET Beam)
  add_subdirectory(Beam)
endif()
if(NOT TARGET Elasticity)
  add_subdirectory(Elasticity)
endif()
if(NOT TARGET LinearElasticity)
  add_subdirectory(Linear)
endif()
if(NOT TARGET Shell)
  add_subdirectory(Shell)
endif()

ifem_add_dependency_dir(
  APP
    Elasticity
  DEPTH
    2
  DIR
    BeamGeometry
    IFEM-BeamGeometry
  TARGET
    IFEM::BeamGeometry
  OPTIONAL
)

ifem_add_application(
  NAME
    LinEl
  SOURCES
    main_LinEl.C
  LIBRARIES
    LinearElasticity
)

if(TARGET IFEM::BeamGeometry)
  target_link_libraries(LinEl PRIVATE IFEM::BeamGeometry)
endif()

ifem_add_application(
  NAME
    ShellEl
  SOURCES
    main_ShellEl.C
  LIBRARIES
    Shell
)

# Regression tests
enable_testing()

ifem_add_hdf5_test(
  TARGET
    LinEl
  TEST_FILES
    Linear/MidShipFrame.hreg
    Linear/ScordelisLo-p4-4x4.hreg
)

ifem_add_hdf5_test(
  TARGET
    LinEl
  CONDITION
    LRSpline_FOUND
  TEST_FILES
    Linear/LR/Lshape-adap-p2.hreg
)

ifem_add_vtf_test(
  TARGET
    LinEl
  TEST_FILES
    Linear/Annulus3D.vreg
    Linear/Hole2D-p3.vreg
    Linear/MidShipFrame.vreg
    Linear/ScordelisLo-p4-4x4.vreg
    Linear/SSUprofil-p1-lagrange.vreg
    Linear/SSUprofil-p1.vreg
    Linear/triangular-neumann.vreg
    Linear/Utkrager.vreg
)

ifem_add_vtf_test(
  TARGET
    LinEl
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    Linear/LR/NavierIcrack_p2-adap.vreg
    Linear/LR/PointLoadedPlate-adap.vreg
)

ifem_add_vtf_test(
  TARGET
    ShellEl
  TEST_FILES
    Shell/PinchedHemisphere-p4.vreg
)

ifem_add_restart_test(
  TARGET
    LinEl
  TEST_FILES
    Linear/CanTS2D-p2-restart.reg
  RESTART_LEVEL
    15
)

ifem_add_restart_test(
  TARGET
    ShellEl
  TEST_FILES
    Shell/Cantilever-p3-restart.reg
  RESTART_LEVEL
    5
)

ifem_add_regression_test(
  TARGET
    LinEl
  TEST_FILES
    Linear/Annulus2D.reg
    Linear/Annulus3D.reg
    Linear/AnnulusWithTemp2D.reg
    Linear/BeamFrame.reg
    Linear/Beam+KLcyl-p3.reg
    Linear/Beam+supel.reg
    Linear/Beam+Uprofil-p2.reg
    Linear/Biharmonic1D-p3.reg
    Linear/boatbeam.reg
    Linear/Bolt.reg
    Linear/CanTD2D-p1.reg
    Linear/Cantilever-Cable3p.reg
    Linear/Cantilever-KLplate.reg
    Linear/Cantilever-KLshell.reg
    Linear/CanTS2D-p1.reg
    Linear/CanTS2D-p2-dmp.reg
    Linear/CanTS2D-p2-dynamic.reg
    Linear/CanTS2D-p2-dyn.reg
    Linear/CanTS2D-p2-qstat.reg
    Linear/CanTS2D-p2.reg
    Linear/CanTS2D-p3.reg
    Linear/Cylinder-Axisymm.reg
    Linear/Cylinder-Lagrange.reg
    Linear/Cylinder-NURBS.reg
    Linear/Cylinder-Spectral.reg
    Linear/EBbeam+KLplate-p2.reg
    Linear/exact_p1.reg
    Linear/exact_p2.reg
    Linear/Harmonic1D-p3.reg
    Linear/Hole2D-immS1-p1.reg
    Linear/Hole2D-immS2-p1.reg
    Linear/Hole2D-Lagrange.reg
    Linear/Hole2D-NURBS.reg
    Linear/Hole2D-Spectral.reg
    Linear/L1-50x50-m.reg
    Linear/MaterialExpression.reg
    Linear/MaterialImage.reg
    Linear/MaterialMultiple.reg
    Linear/MidShipFrame.reg
    Linear/NavierPart_p2.reg
    Linear/NavierPart_p3.reg
    Linear/NavierPoint_p3.reg
    Linear/NavierUload_p2.reg
    Linear/NavierUload_p3.reg
    Linear/NREL_5MW_blade-vibration.reg
    Linear/PinchedHemisphere-p4.reg
    Linear/PipeJoint-Lagrange.reg
    Linear/PipeJoint-NURBS.reg
    Linear/PipeJoint-vibration.reg
    Linear/PlateLineLoad.reg
    Linear/PlateWithColumns_p2.reg
    Linear/PointLoadedPlate-project-file.reg
    Linear/PointLoadedPlate-project.reg
    Linear/PolePinchedSphere-p3.reg
    Linear/PressurizedHemisphere-p2.reg
    Linear/RectPlate-dynamic.reg
    Linear/RectPlate-modal.reg
    Linear/ScordelisLo-p4-4x4.reg
    Linear/ScordelisLo-symm-p4-2x2.reg
    Linear/SquarePlate-Lagrange.reg
    Linear/SquarePlate-Spectral.reg
    Linear/SquarePlate-Splines.reg
    Linear/SquarePoint_p3.reg
    Linear/SquarePoint-varying.reg
    Linear/SS45-Cable4p.reg
    Linear/SSbeamPointLoad-adap.reg
    Linear/SSbeamPointLoad-graded.reg
    Linear/SSbeamPointLoad-project-file.reg
    Linear/SSbeamPointLoad-project.reg
    Linear/SSbeamPointLoad.reg
    Linear/SScablePointLoad.reg
    Linear/SSmembrane-p1.reg
    Linear/SSsolid-p2.reg
    Linear/SSUprofil-p1.reg
    Linear/stretch-p1.reg
    Linear/triangular-neumann.reg
    Linear/Tripod.reg
    Linear/Utkrager.reg
    Linear/Utkrager-Box.reg
    Linear/Utkrager-lumped.reg
    Linear/Utkrager-dyn.reg
    Linear/Utkrager-Ecc.reg
    Linear/Utkrager-Pipe.reg
    Linear/Utkrager-dump.reg
    Linear/Utkrager-load.reg
    Linear/WedgeStretch.reg
)

ifem_add_regression_test(
  TARGET
    ShellEl
  TEST_FILES
    Shell/Cantilever-p2-1x1.reg
    Shell/Cantilever-p3-10x1.reg
    Shell/HingedShallowArch-p2.reg
    Shell/HingedShallowArch-press.reg
    Shell/PinchedHemisphere-p4.reg
    Shell/PolePinchedSphere_p2.reg
)

if(NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
  set_tests_properties(LinEl+Linear/Utkrager-load.reg
                       PROPERTIES DEPENDS
                       LinEl+Linear/Utkrager-dump.reg)

  set_tests_properties(LinEl+Linear/Utkrager-dump.reg
                       PROPERTIES FIXTURES_SETUP Dump)
  set_tests_properties(LinEl+Linear/Utkrager-load.reg
                       PROPERTIES FIXTURES_REQUIRED Dump)
endif()

ifem_add_regression_test(
  TARGET
    LinEl
  CONDITIONS
    HDF5_FOUND
  TEST_FILES
    Linear/CanTS2D-fields-p2.reg
)

ifem_add_regression_test(
  TARGET
    LinEl
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    Linear/LR/Cantilever-KLshell.reg
    Linear/LR/CanTS2D-p1.reg
    Linear/LR/CanTS2D-p1-VCPs.reg
    Linear/LR/disc-stair-adap.reg
    Linear/LR/exact_p2.reg
    Linear/LR/Lshape-adap-p2.reg
    Linear/LR/NavierIcrack_p2-adap.reg
    Linear/LR/NavierIcrack_p2.reg
    Linear/LR/PointLoadedPlate-adap.reg
    Linear/LR/PointLoadedPlate-NURBS.reg
    Linear/LR/PointLoadedPlate+project.reg
    Linear/LR/ScordelisPoint-NURBS.reg
    Linear/LR/stretch-p1.reg
)

ifem_add_regression_test(
  TARGET
    LinEl
  CONDITIONS
    SPR_FOUND
  TEST_FILES
    Linear/SPR/Annulus2D.reg
    Linear/SPR/Annulus3D.reg
    Linear/SPR/AnnulusWithTemp2D.reg
    Linear/SPR/BeamFrame.reg
    Linear/SPR/Beam+KLcyl-p3.reg
    Linear/SPR/Beam+supel.reg
    Linear/SPR/Beam+Uprofil-p2.reg
    Linear/SPR/Biharmonic1D-p3.reg
    Linear/SPR/boatbeam.reg
    Linear/SPR/CanTD2D-p1.reg
    Linear/SPR/Cantilever-Cable3p.reg
    Linear/SPR/Cantilever-KLplate.reg
    Linear/SPR/Cantilever-KLshell.reg
    Linear/SPR/CanTS2D-p1.reg
    Linear/SPR/CanTS2D-p2-dmp.reg
    Linear/SPR/CanTS2D-p2-dynamic.reg
    Linear/SPR/CanTS2D-p2-dyn.reg
    Linear/SPR/CanTS2D-p2-qstat.reg
    Linear/SPR/CanTS2D-p2.reg
    Linear/SPR/CanTS2D-p3.reg
    Linear/SPR/Cylinder-Axisymm.reg
    Linear/SPR/Cylinder-Lagrange.reg
    Linear/SPR/Cylinder-NURBS.reg
    Linear/SPR/Cylinder-Spectral.reg
    Linear/SPR/EBbeam+KLplate-p2.reg
    Linear/SPR/exact_p1.reg
    Linear/SPR/exact_p2.reg
    Linear/SPR/Harmonic1D-p3.reg
    Linear/SPR/Hole2D-immS1-p1.reg
    Linear/SPR/Hole2D-immS2-p1.reg
    Linear/SPR/Hole2D-Lagrange.reg
    Linear/SPR/Hole2D-NURBS.reg
    Linear/SPR/Hole2D-Spectral.reg
    Linear/SPR/L1-50x50-m.reg
    Linear/SPR/MaterialExpression.reg
    Linear/SPR/MaterialImage.reg
    Linear/SPR/MaterialMultiple.reg
    Linear/SPR/MidShipFrame.reg
    Linear/SPR/NavierPart_p2.reg
    Linear/SPR/NavierPart_p3.reg
    Linear/SPR/NavierPoint_p3.reg
    Linear/SPR/NavierUload_p2.reg
    Linear/SPR/NavierUload_p3.reg
    Linear/SPR/NREL_5MW_blade-vibration.reg
    Linear/SPR/PinchedHemisphere-p4.reg
    Linear/SPR/PipeJoint-Lagrange.reg
    Linear/SPR/PipeJoint-NURBS.reg
    Linear/SPR/PlateLineLoad.reg
    Linear/SPR/PlateWithColumns_p2.reg
    Linear/SPR/PointLoadedPlate-project-file.reg
    Linear/SPR/PointLoadedPlate-project.reg
    Linear/SPR/PolePinchedSphere-p3.reg
    Linear/SPR/RectPlate-dynamic.reg
    Linear/SPR/RectPlate-modal.reg
    Linear/SPR/ScordelisLo-p4-4x4.reg
    Linear/SPR/ScordelisLo-symm-p4-2x2.reg
    Linear/SPR/SquarePlate-Lagrange.reg
    Linear/SPR/SquarePlate-Spectral.reg
    Linear/SPR/SquarePlate-Splines.reg
    Linear/SPR/SquarePoint_p3.reg
    Linear/SPR/SS45-Cable4p.reg
    Linear/SPR/SSbeamPointLoad-adap.reg
    Linear/SPR/SSbeamPointLoad-graded.reg
    Linear/SPR/SSbeamPointLoad-project-file.reg
    Linear/SPR/SSbeamPointLoad-project.reg
    Linear/SPR/SSbeamPointLoad.reg
    Linear/SPR/SScablePointLoad.reg
    Linear/SPR/SSmembrane-p1.reg
    Linear/SPR/SSsolid-p2.reg
    Linear/SPR/SSUprofil-p1.reg
    Linear/SPR/stretch-p1.reg
    Linear/SPR/triangular-neumann.reg
    Linear/SPR/Utkrager.reg
    Linear/SPR/WedgeStretch.reg
)

# Unit tests
ifem_add_test_app(
  NAME
    LinEl
  SOURCES
    Linear/Test/TestStaticCondensation.C
  WORKDIR
    ${PROJECT_SOURCE_DIR}/Test/Linear
  LIBRARIES
    LinearElasticity
)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  ifem_add_check_target()
endif()

# For generating the doxy
ifem_add_doc_target(
  TARGET
    Elasticity
  EXTRA_DIRS
     ${BeamGeometry_DIRECTORY}
)

# Installation
include(GNUInstallDirs)
install(TARGETS LinEl ShellEl DESTINATION ${CMAKE_INSTALL_BINDIR})
