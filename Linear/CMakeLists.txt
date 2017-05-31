project(LinEl)

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

include_directories(${IFEM_INCLUDES} ${PROJECT_SOURCE_DIR}/..
                                     ${PROJECT_SOURCE_DIR}/../Beam)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

if(NOT TARGET Elasticity)
  add_subdirectory(../ Elasticity)
endif()
if(NOT TARGET Beam)
  add_subdirectory(../Beam Beam)
endif()

file(GLOB LinEl_SRCS *.C)

add_executable(LinEl ${LinEl_SRCS})
target_link_libraries(LinEl Beam Elasticity ${IFEM_LIBRARIES})

list(APPEND CHECK_SOURCES ${LinEl_SRCS})

# Installation
install(TARGETS LinEl DESTINATION bin COMPONENT bin)

# For generating the doxy
set(EXTRA_DOXY_PATHS "${PROJECT_SOURCE_DIR} \\
                      ${PROJECT_BINARY_DIR} \\
                      ${PROJECT_SOURCE_DIR}/.. \\
                      ${PROJECT_SOURCE_DIR}/../Beam")
add_doc_target(LinEl LinEl)

if(IFEM_INSTALL_DOXY)
  install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_BUILD_TOOL} doc WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\")" COMPONENT doc)
  install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html DESTINATION ${CMAKE_INSTALL_DOCDIR}/Apps/LinEl
          COMPONENT doc
          PATTERN *.md5 EXCLUDE
          PATTERN *.map EXCLUDE)
endif(IFEM_INSTALL_DOXY)

# Regression tests
enable_testing()
include(IFEMTesting)

if(MPI_FOUND)
  # Add parallel tests here
else()
  file(GLOB LINEL_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                                     ${PROJECT_SOURCE_DIR}/Test/*.reg)
  if(IFEM_USE_LRSPLINES AND LRSpline_FOUND)
    file(GLOB LRE_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                                     ${PROJECT_SOURCE_DIR}/Test/LR/*.reg)
  endif()
  foreach(TESTFILE ${LINEL_TESTFILES} ${LRE_TESTFILES})
    ifem_add_test(${TESTFILE} LinEl)
  endforeach()
endif()

list(APPEND TEST_APPS LinEl)
if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  # Remove when there are unit tests in this application
  if(NOT UNIT_TEST_NUMBER)
    set(UNIT_TEST_NUMBER 0)
  endif()
  add_check_target()
endif()
