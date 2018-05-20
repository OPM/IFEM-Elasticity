project(ShellEl)

cmake_minimum_required(VERSION 2.6)

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
  foreach(TESTFILE ${TESTFILES})
    ifem_add_test(${TESTFILE} ShellEl)
  endforeach()
endif()

list(APPEND TEST_APPS ShellEl)
if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else(IFEM_COMMON_APP_BUILD)
  # Remove when there are unit tests in this application
  if(NOT UNIT_TEST_NUMBER)
    set(UNIT_TEST_NUMBER 0)
  endif()
  add_check_target()
endif()
