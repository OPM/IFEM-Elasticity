cmake_minimum_required(VERSION 3.5)

project(ShellEl)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../../cmake/Modules)

# Required packages
if(NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
  include_directories(${IFEM_INCLUDE_DIRS})
  add_definitions(${IFEM_DEFINITIONS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")
endif()

include_directories(${IFEM_INCLUDES} ${PROJECT_SOURCE_DIR}/..)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

if(NOT TARGET Elasticity)
  add_subdirectory(.. Elasticity)
endif()

file(GLOB Shell_SRCS *.C)

add_executable(ShellEl ${Shell_SRCS})
target_link_libraries(ShellEl Elasticity ${IFEM_LIBRARIES})

# Installation

install(TARGETS ShellEl DESTINATION bin COMPONENT bin)
list(APPEND CHECK_SOURCES ${Shell_SRCS})

# For generating the doxy
set(EXTRA_DOXY_PATHS "${PROJECT_SOURCE_DIR} \\
                      ${PROJECT_SOURCE_DIR}/..")
add_doc_target(ShellEl ShellEl)

# Regression tests
enable_testing()
include(IFEMTesting)

if(NOT MPI_FOUND)
  # Only serial tests for this app
  file(GLOB TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                      ${PROJECT_SOURCE_DIR}/Test/*.reg)
  list(REMOVE_ITEM TESTFILES Cantilever-p3-restart.reg)
  file(GLOB VTF_FILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                      ${PROJECT_SOURCE_DIR}/Test/*.vreg)
  foreach(TESTFILE ${TESTFILES})
    ifem_add_test(${TESTFILE} ShellEl)
  endforeach()
  foreach(VTF_FILE ${VTF_FILES})
    ifem_add_vtf_test(${VTF_FILE} ShellEl)
  endforeach()
  if(HDF5_FOUND AND CEREAL_FOUND)
    ifem_add_restart_test(Cantilever-p3-restart.reg ShellEl 5)
  endif()
endif()

list(APPEND TEST_APPS ShellEl)
if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else(IFEM_COMMON_APP_BUILD)
  add_check_target()
endif()
