cmake_minimum_required(VERSION 3.22)

project(FiniteDeformation LANGUAGES C CXX Fortran)

ifem_add_library(
  NAME
    FiniteDeformation
  ALIAS
    IFEM::FiniteDeformation
  SOURCES
    DruckerPrager.C
    LinearMaterial.C
    MixedTanMat.C
    MortarContact.C
    NeoHookeMaterial.C
    NLoptions.C
    NonlinearElasticity.C
    NonlinearElasticityFbar.C
    NonlinearElasticityTL.C
    NonlinearElasticityUL.C
    NonlinearElasticityULMixed.C
    NonlinearElasticityULMX.C
    PlasticMaterial.C
    RigidBody.C
    SIMContact.C
    SIMFiniteDefEl.C
    accKM.f
    accKMx.f
    stiff_TL.f
  HEADERS
    DruckerPrager.h
    LinearMaterial.h
    MixedTanMat.h
    MortarContact.h
    NeoHookeMaterial.h
    NLoptions.h
    NonlinearElasticityFbar.h
    NonlinearElasticity.h
    NonlinearElasticityTL.h
    NonlinearElasticityUL.h
    NonlinearElasticityULMixed.h
    NonlinearElasticityULMX.h
    PlasticMaterial.h
    RigidBody.h
    SIMContact.h
    SIMFiniteDefEl.h
    SIMLinEl.h
  LIBRARIES
    IFEM::Elasticity
)

target_compile_definitions(FiniteDeformation PUBLIC USE_FTNMAT=1)
