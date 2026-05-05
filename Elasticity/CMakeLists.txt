cmake_minimum_required(VERSION 3.22)

project(Elasticity LANGUAGES C CXX)

ifem_add_library(
  NAME
    Elasticity
  SOURCES
    ArcLengthDriver.C
    ElasticBase.C
    ElasticityArgs.C
    Elasticity.C
    IsotropicTextureMat.C
    KirchhoffLove.C
    LinearElasticity.C
    LinIsotropic.C
    LocalSystems.C
    NonlinearDriver.C
    SIMElasticity.C
    SIMElasticityWrap.C
    SIMRigid.C
  HEADERS
    ArcLengthDriver.h
    ElasticBase.h
    ElasticityArgs.h
    Elasticity.h
    ElasticityUtils.h
    IsotropicTextureMat.h
    KirchhoffLove.h
    LinearElasticity.h
    LinIsotropic.h
    MaterialBase.h
    NewmarkDriver.h
    NonlinearDriver.h
    SIMElasticity.h
    SIMElasticityWrap.h
    SIMRigid.h
  LIBRARIES
    IFEM
)
