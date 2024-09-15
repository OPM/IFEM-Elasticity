cmake_minimum_required(VERSION 3.5)

project(LinEl)
message(STATUS "Configuring project ${PROJECT_NAME}")

if(NOT IFEM_CONFIGURED)
  list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../../cmake/Modules)
  find_package(IFEM REQUIRED)
  include_directories(${IFEM_INCLUDE_DIRS})
  add_definitions(${IFEM_DEFINITIONS})
  string(APPEND CMAKE_CXX_FLAGS " ${IFEM_CXX_FLAGS}")
endif()

include_directories(${IFEM_INCLUDES} ${PROJECT_SOURCE_DIR}
                                     ${PROJECT_SOURCE_DIR}/..
                                     ${PROJECT_SOURCE_DIR}/../Shell
                                     ${PROJECT_SOURCE_DIR}/../Beam)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

if(NOT TARGET Elasticity)
  add_subdirectory(../ Elasticity)
endif()
if(NOT TARGET Beam)
  add_subdirectory(../Beam Beam)
endif()
if(EXISTS ${PROJECT_SOURCE_DIR}/../IFEM-BeamGeometry)
  set(BEAMGEO_SRC ../IFEM-BeamGeometry)
else()
  set(BEAMGEO_SRC ../BeamGeometry)
endif()
if(NOT TARGET BeamGeometry AND EXISTS ${PROJECT_SOURCE_DIR}/${BEAMGEO_SRC})
  message(STATUS "Including ${PROJECT_SOURCE_DIR}/${BEAMGEO_SRC} for project ${CMAKE_PROJECT_NAME}")
  add_subdirectory(${BEAMGEO_SRC} BeamGeometry)
  include_directories(${PROJECT_SOURCE_DIR}/${BEAMGEO_SRC})
  string(APPEND CMAKE_CXX_FLAGS " -DHAS_GEOMETRY")
endif()

file(GLOB LinEl_SRCS *.C)
list(APPEND LinEl_SRCS ../Shell/SIMKLShell.C ../Shell/KirchhoffLoveShell.C)

add_executable(LinEl ${LinEl_SRCS})
target_link_libraries(LinEl Beam Elasticity ${IFEM_LIBRARIES})
if(TARGET BeamGeometry)
  target_link_libraries(LinEl BeamGeometry)
endif()

list(APPEND CHECK_SOURCES ${LinEl_SRCS})

# Installation
install(TARGETS LinEl DESTINATION bin COMPONENT bin)

# For generating the doxy
set(EXTRA_DOXY_PATHS "${PROJECT_SOURCE_DIR} \\
                      ${PROJECT_BINARY_DIR} \\
                      ${PROJECT_SOURCE_DIR}/.. \\
                      ${PROJECT_SOURCE_DIR}/../Beam")
if(EXISTS ${PROJECT_SOURCE_DIR}/${BEAMGEO_SRC})
  string(APPEND EXTRA_DOXY_PATHS " ${PROJECT_SOURCE_DIR}/${BEAMGEO_SRC}")
endif()
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

file(GLOB LINEL_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                                   ${PROJECT_SOURCE_DIR}/Test/*.reg)
file(GLOB LINEL_VTF_FILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                                   ${PROJECT_SOURCE_DIR}/Test/*.vreg)
file(GLOB LINEL_HDF5FILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                                   ${PROJECT_SOURCE_DIR}/Test/*.hreg)
if(LRSpline_FOUND OR LRSPLINE_FOUND)
  file(GLOB LRE_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                                   ${PROJECT_SOURCE_DIR}/Test/LR/*.reg)
  file(GLOB LRE_VTF_FILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                                   ${PROJECT_SOURCE_DIR}/Test/LR/*.vreg)
  file(GLOB LRE_HDF5FILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                                   ${PROJECT_SOURCE_DIR}/Test/LR/*.hreg)
endif()
if(SPR_LIBRARIES)
  file(GLOB SPR_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                                   ${PROJECT_SOURCE_DIR}/Test/SPR/*.reg)
endif()

list(REMOVE_ITEM LINEL_TESTFILES CanTS2D-p2-restart.reg)

if(NOT HDF5_FOUND)
  list(REMOVE_ITEM LINEL_TESTFILES CanTS2D-fields-p2.reg)
endif()

if(NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
  foreach(TESTFILE ${LINEL_TESTFILES} ${LRE_TESTFILES} ${SPR_TESTFILES})
    ifem_add_test(${TESTFILE} LinEl)
  endforeach()
  foreach(TESTFILE ${LINEL_VTF_FILES} ${LRE_VTF_FILES})
    ifem_add_vtf_test(${TESTFILE} LinEl)
  endforeach()
  foreach(TESTFILE ${LINEL_HDF5FILES} ${LRE_HDF5FILES})
    ifem_add_hdf5_test(${TESTFILE} LinEl)
  endforeach()
  if(HDF5_FOUND AND CEREAL_FOUND)
    ifem_add_restart_test(CanTS2D-p2-restart.reg LinEl 15)
  endif()
endif()

list(APPEND TEST_APPS LinEl)

# Unit tests
add_library(LinearElasticity STATIC ${PROJECT_SOURCE_DIR}/SIMLinEl2D.C ${PROJECT_SOURCE_DIR}/AnalyticSolutions.C)
IFEM_add_test_app(${PROJECT_SOURCE_DIR}/Test/*.C
                  ${PROJECT_SOURCE_DIR}/Test
                  LinEl 0
                  LinearElasticity Elasticity ${IFEM_LIBRARIES})
list(APPEND TEST_APPS LinEl-test)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else()
  add_check_target()
endif()
