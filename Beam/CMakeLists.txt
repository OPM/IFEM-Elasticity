cmake_minimum_required(VERSION 3.5)

project(Beam)

include_directories(${PROJECT_SOURCE_DIR}/..)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

file(GLOB Beam_SOURCES *.C)
file(GLOB Beam_HEADERS *.h)
add_library(Beam STATIC ${Beam_SOURCES})

list(APPEND CHECK_SOURCES ${Beam_SOURCES})
set(CHECK_SOURCES ${CHECK_SOURCES} PARENT_SCOPE)
