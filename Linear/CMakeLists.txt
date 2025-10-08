cmake_minimum_required(VERSION 3.22)

project(LinEl LANGUAGES C CXX)

message(STATUS "Configuring project ${PROJECT_NAME}")

if(NOT TARGET IFEM)
  if(IFEM_AS_SUBMODULE)
    add_subdirectory(../../.. IFEM)
    set(IFEM_PATH ${PROJECT_SOURCE_DIR}/../../..)
  else()
    get_filename_component(BUILD_DIR ${PROJECT_BINARY_DIR} NAME)
    set(IFEM_DIR ${PROJECT_SOURCE_DIR}/../../../${BUILD_DIR})
    find_package(IFEM REQUIRED)
  endif()
endif()

if(NOT TARGET Elasticity)
  add_subdirectory(../ Elasticity)
endif()
if(NOT TARGET Beam)
  add_subdirectory(../Beam Beam)
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
    BeamGeometry
  OPTIONAL
)

ifem_add_library(
  NAME
    LinearElasticity
  BASE_DIRS
    ${PROJECT_SOURCE_DIR}/../Shell
  SOURCES
    AnalyticSolutions.C
    DynamicDriver.C
    KirchhoffLovePlate.C
    SIMLinEl2D.C
    ModalDriver.C
    MultiLoadCaseDriver.C
    SIMLinEl2D.C
    SIMLinEl3D.C
    SIMLinElBeamC1.C
    SIMLinElKL.C
    SIMLinElSup.C
    SIMLinKLModal.C
    SIMmcStatic.C
    ../Shell/SIMKLShell.C
    ../Shell/KirchhoffLoveShell.C
  HEADERS
    AnalyticSolutions.h
    DynamicSim.h
    KirchhoffLovePlate.h
    ModalDriver.h
    MultiLoadCaseSim.h
    SIMLinElBeamC1.h
    SIMLinEl.h
    SIMLinElKL.h
    SIMLinElModal.h
    SIMLinElSup.h
    SIMLinKLModal.h
    SIMmcStatic.h
    ../Shell/SIMKLShell.h
    ../Shell/KirchhoffLoveShell.h
  LIBRARIES
    Beam
    Elasticity
)

ifem_add_application(
  NAME
    LinEl
  SOURCES
    main_LinEl.C
  LIBRARIES
    LinearElasticity
)

if(TARGET BeamGeometry)
  target_link_libraries(LinEl PRIVATE BeamGeometry)
endif()

# Regression tests
enable_testing()

ifem_add_hdf5_test(
  TARGET
    LinEl
  TEST_FILES
    MidShipFrame.hreg
    ScordelisLo-p4-4x4.hreg
)

ifem_add_hdf5_test(
  TARGET
    LinEl
  CONDITION
    LRSpline_FOUND
  TEST_FILES
    LR/Lshape-adap-p2.hreg
)

ifem_add_vtf_test(
  TARGET
    LinEl
  TEST_FILES
    Annulus3D.vreg
    Hole2D-p3.vreg
    MidShipFrame.vreg
    ScordelisLo-p4-4x4.vreg
    SSUprofil-p1-lagrange.vreg
    SSUprofil-p1.vreg
    triangular-neumann.vreg
    Utkrager.vreg
)

ifem_add_vtf_test(
  TARGET
    LinEl
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    LR/NavierIcrack_p2-adap.vreg
    LR/PointLoadedPlate-adap.vreg
)

ifem_add_restart_test(
  TARGET
    LinEl
  TEST_FILES
    CanTS2D-p2-restart.reg
  RESTART_LEVEL
    15
)

ifem_add_regression_test(
  TARGET
    LinEl
  TEST_FILES
    Annulus2D.reg
    Annulus3D.reg
    AnnulusWithTemp2D.reg
    BeamFrame.reg
    Beam+KLcyl-p3.reg
    Beam+supel.reg
    Beam+Uprofil-p2.reg
    Biharmonic1D-p3.reg
    boatbeam.reg
    Bolt.reg
    CanTD2D-p1.reg
    Cantilever-Cable3p.reg
    Cantilever-KLplate.reg
    Cantilever-KLshell.reg
    CanTS2D-p1.reg
    CanTS2D-p2-dmp.reg
    CanTS2D-p2-dynamic.reg
    CanTS2D-p2-dyn.reg
    CanTS2D-p2-qstat.reg
    CanTS2D-p2.reg
    CanTS2D-p3.reg
    Cylinder-Axisymm.reg
    Cylinder-Lagrange.reg
    Cylinder-NURBS.reg
    Cylinder-Spectral.reg
    EBbeam+KLplate-p2.reg
    exact_p1.reg
    exact_p2.reg
    Harmonic1D-p3.reg
    Hole2D-immS1-p1.reg
    Hole2D-immS2-p1.reg
    Hole2D-Lagrange.reg
    Hole2D-NURBS.reg
    Hole2D-Spectral.reg
    L1-50x50-m.reg
    MaterialExpression.reg
    MaterialImage.reg
    MaterialMultiple.reg
    MidShipFrame.reg
    NavierPart_p2.reg
    NavierPart_p3.reg
    NavierPoint_p3.reg
    NavierUload_p2.reg
    NavierUload_p3.reg
    NREL_5MW_blade-vibration.reg
    PinchedHemisphere-p4.reg
    PipeJoint-Lagrange.reg
    PipeJoint-NURBS.reg
    PipeJoint-vibration.reg
    PlateLineLoad.reg
    PlateWithColumns_p2.reg
    PointLoadedPlate-project-file.reg
    PointLoadedPlate-project.reg
    PolePinchedSphere-p3.reg
    PressurizedHemisphere-p2.reg
    RectPlate-dynamic.reg
    RectPlate-modal.reg
    ScordelisLo-p4-4x4.reg
    ScordelisLo-symm-p4-2x2.reg
    SquarePlate-Lagrange.reg
    SquarePlate-Spectral.reg
    SquarePlate-Splines.reg
    SquarePoint_p3.reg
    SquarePoint-varying.reg
    SS45-Cable4p.reg
    SSbeamPointLoad-adap.reg
    SSbeamPointLoad-graded.reg
    SSbeamPointLoad-project-file.reg
    SSbeamPointLoad-project.reg
    SSbeamPointLoad.reg
    SScablePointLoad.reg
    SSmembrane-p1.reg
    SSsolid-p2.reg
    SSUprofil-p1.reg
    stretch-p1.reg
    triangular-neumann.reg
    Tripod.reg
    Utkrager.reg
    Utkrager-Box.reg
    Utkrager-lumped.reg
    Utkrager-dyn.reg
    Utkrager-Ecc.reg
    Utkrager-Pipe.reg
    Utkrager-dump.reg
    Utkrager-load.reg
    WedgeStretch.reg
)

if(NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
  set_tests_properties(LinEl+Utkrager-load.reg
                       PROPERTIES DEPENDS
                       LinEl+Utkrager-dump.reg)

  set_tests_properties(LinEl+Utkrager-dump.reg
                       PROPERTIES FIXTURES_SETUP Dump)
  set_tests_properties(LinEl+Utkrager-load.reg
                       PROPERTIES FIXTURES_REQUIRED Dump)
endif()

ifem_add_regression_test(
  TARGET
    LinEl
  CONDITIONS
    HDF5_FOUND
  TEST_FILES
    CanTS2D-fields-p2.reg
)

ifem_add_regression_test(
  TARGET
    LinEl
  CONDITIONS
    LRSpline_FOUND
  TEST_FILES
    LR/Cantilever-KLshell.reg
    LR/CanTS2D-p1.reg
    LR/CanTS2D-p1-VCPs.reg
    LR/disc-stair-adap.reg
    LR/exact_p2.reg
    LR/Lshape-adap-p2.reg
    LR/NavierIcrack_p2-adap.reg
    LR/NavierIcrack_p2.reg
    LR/PointLoadedPlate-adap.reg
    LR/PointLoadedPlate-NURBS.reg
    LR/PointLoadedPlate+project.reg
    LR/ScordelisPoint-NURBS.reg
    LR/stretch-p1.reg
)

ifem_add_regression_test(
  TARGET
    LinEl
  CONDITIONS
    SPR_FOUND
  TEST_FILES
    SPR/Annulus2D.reg
    SPR/Annulus3D.reg
    SPR/AnnulusWithTemp2D.reg
    SPR/BeamFrame.reg
    SPR/Beam+KLcyl-p3.reg
    SPR/Beam+supel.reg
    SPR/Beam+Uprofil-p2.reg
    SPR/Biharmonic1D-p3.reg
    SPR/boatbeam.reg
    SPR/CanTD2D-p1.reg
    SPR/Cantilever-Cable3p.reg
    SPR/Cantilever-KLplate.reg
    SPR/Cantilever-KLshell.reg
    SPR/CanTS2D-p1.reg
    SPR/CanTS2D-p2-dmp.reg
    SPR/CanTS2D-p2-dynamic.reg
    SPR/CanTS2D-p2-dyn.reg
    SPR/CanTS2D-p2-qstat.reg
    SPR/CanTS2D-p2.reg
    SPR/CanTS2D-p3.reg
    SPR/Cylinder-Axisymm.reg
    SPR/Cylinder-Lagrange.reg
    SPR/Cylinder-NURBS.reg
    SPR/Cylinder-Spectral.reg
    SPR/EBbeam+KLplate-p2.reg
    SPR/exact_p1.reg
    SPR/exact_p2.reg
    SPR/Harmonic1D-p3.reg
    SPR/Hole2D-immS1-p1.reg
    SPR/Hole2D-immS2-p1.reg
    SPR/Hole2D-Lagrange.reg
    SPR/Hole2D-NURBS.reg
    SPR/Hole2D-Spectral.reg
    SPR/L1-50x50-m.reg
    SPR/MaterialExpression.reg
    SPR/MaterialImage.reg
    SPR/MaterialMultiple.reg
    SPR/MidShipFrame.reg
    SPR/NavierPart_p2.reg
    SPR/NavierPart_p3.reg
    SPR/NavierPoint_p3.reg
    SPR/NavierUload_p2.reg
    SPR/NavierUload_p3.reg
    SPR/NREL_5MW_blade-vibration.reg
    SPR/PinchedHemisphere-p4.reg
    SPR/PipeJoint-Lagrange.reg
    SPR/PipeJoint-NURBS.reg
    SPR/PlateLineLoad.reg
    SPR/PlateWithColumns_p2.reg
    SPR/PointLoadedPlate-project-file.reg
    SPR/PointLoadedPlate-project.reg
    SPR/PolePinchedSphere-p3.reg
    SPR/RectPlate-dynamic.reg
    SPR/RectPlate-modal.reg
    SPR/ScordelisLo-p4-4x4.reg
    SPR/ScordelisLo-symm-p4-2x2.reg
    SPR/SquarePlate-Lagrange.reg
    SPR/SquarePlate-Spectral.reg
    SPR/SquarePlate-Splines.reg
    SPR/SquarePoint_p3.reg
    SPR/SS45-Cable4p.reg
    SPR/SSbeamPointLoad-adap.reg
    SPR/SSbeamPointLoad-graded.reg
    SPR/SSbeamPointLoad-project-file.reg
    SPR/SSbeamPointLoad-project.reg
    SPR/SSbeamPointLoad.reg
    SPR/SScablePointLoad.reg
    SPR/SSmembrane-p1.reg
    SPR/SSsolid-p2.reg
    SPR/SSUprofil-p1.reg
    SPR/stretch-p1.reg
    SPR/triangular-neumann.reg
    SPR/Utkrager.reg
    SPR/WedgeStretch.reg
)

# Unit tests
ifem_add_test_app(
  NAME
    LinEl
  SOURCES
    Test/TestStaticCondensation.C
  WORKDIR
    ${PROJECT_SOURCE_DIR}/Test
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
    LinEl
  EXTRA_DIRS
     ${BeamGeometry_DIRECTORY}
)

# Installation
include(GNUInstallDirs)
install(TARGETS LinEl DESTINATION ${CMAKE_INSTALL_BINDIR})
