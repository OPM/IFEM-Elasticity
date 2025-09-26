cmake_minimum_required(VERSION 3.22)

project(Beam LANGUAGES C CXX)

ifem_add_library(
  NAME
    Beam
  SOURCES
    BeamProperty.C
    ElasticBar.C
    ElasticBeam.C
    ElasticCable.C
    SIMElastic1D.C
    SIMElasticBar.C
  HEADERS
    BeamProperty.h
    ElasticBar.h
    ElasticBeam.h
    ElasticCable.h
    SIMElastic1D.h
    SIMElasticBar.h
    SIMLinElBeam.h
  LIBRARIES
    Elasticity
)

target_include_directories(Beam INTERFACE ${PROJECT_SOURCE_DIR})
