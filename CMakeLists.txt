project(Elasticity)

cmake_minimum_required(VERSION 2.6)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

file(GLOB El_SOURCES *.C)
file(GLOB El_HEADERS *.h)
add_library(Elasticity STATIC ${El_SOURCES})
