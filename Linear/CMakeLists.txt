PROJECT(LinEl)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../../cmake/Modules)

IF (NOT IFEM_CONFIGURED)
  FIND_PACKAGE(IFEM REQUIRED)
ENDIF(NOT IFEM_CONFIGURED)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")

INCLUDE_DIRECTORIES(${IFEM_INCLUDES} ${PROJECT_SOURCE_DIR}/..
                                     ${PROJECT_SOURCE_DIR}/../FiniteDeformation)

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# Generate regtest script with correct paths
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/../../Common/scripts/regtest.sh.in
               regtest.sh)

FILE(GLOB LinEl_SRCS *.C;../*.C;
                     ../FiniteDeformation/SIMElasticBar.C;
                     ../FiniteDeformation/Elastic*.C)
ADD_EXECUTABLE(LinEl ${LinEl_SRCS})
TARGET_LINK_LIBRARIES(LinEl ${IFEM_LIBRARIES})

# Installation
INSTALL(TARGETS LinEl DESTINATION bin COMPONENT bin)

# Regression tests
ENABLE_TESTING()

# Generate regtest script with correct paths
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/../../Common/scripts/regtest.sh.in
               regtest.sh)
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/../../Common/scripts/doregtests.sh.in
               doregtests.sh @ONLY)
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/../../Common/scripts/check-patchseries.sh.in
               check-patchseries.sh)

IF(IFEM_USE_PARALLEL_PETSC)
  # Add parallel tests here
ELSE(IFEM_USE_PARALLEL_PETSC)
  FILE(GLOB_RECURSE LINEL_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/Test
                                    ${PROJECT_SOURCE_DIR}/Test/*.reg)
  FOREACH(TESTFILE ${LINEL_TESTFILES})
    IFEM_ADD_TEST(${TESTFILE} LinEl)
  ENDFOREACH()
ENDIF(IFEM_USE_PARALLEL_PETSC)
