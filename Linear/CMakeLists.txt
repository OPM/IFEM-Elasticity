cmake_minimum_required(VERSION 3.22)

project(LinearElasticity LANGUAGES C CXX)

ifem_add_library(
  NAME
    LinearElasticity
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
  LIBRARIES
    Beam
    Elasticity
    Shell
)
