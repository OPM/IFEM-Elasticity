cmake_minimum_required(VERSION 3.22)

project(ShellEl LANGUAGES C CXX)

if(NOT TARGET IFEM)
  if(IFEM_AS_SUBMODULE)
    add_subdirectory(../../.. IFEM)
    set(IFEM_PATH ${PROJECT_SOURCE_DIR}/../../..)
  else()
    get_filename_component(BUILD_DIR ${PROJECT_BINARY_DIR} NAME)
    set(IFEM_DIR ${PROJECT_SOURCE_DIR}/../../../${BUILD_DIR})
    find_package(IFEM REQUIRED)
  endif()
endif()

if(NOT TARGET Elasticity)
  add_subdirectory(.. Elasticity)
endif()

ifem_add_application(
  NAME
    ShellEl
  SOURCES
    KirchhoffLoveShell.C
    main.C
    NLKirchhoffLoveShell.C
    SIMKLShell.C
    SIMShell.C
  HEADERS
    KirchhoffLoveShell.h
    NLKirchhoffLoveShell.h
    SIMKLShell.h
    SIMShell.h
  LIBRARIES
    Elasticity
)

# Regression tests
enable_testing()

ifem_add_vtf_test(
  TARGET
    ShellEl
  TEST_FILES
    PinchedHemisphere-p4.vreg
)

ifem_add_restart_test(
  TARGET
    ShellEl
  TEST_FILES
    Cantilever-p3-restart.reg
  RESTART_LEVEL
    5
)

ifem_add_regression_test(
  TARGET
    ShellEl
  TEST_FILES
    Cantilever-p2-1x1.reg
    Cantilever-p3-10x1.reg
    HingedShallowArch-p2.reg
    HingedShallowArch-press.reg
    PinchedHemisphere-p4.reg
    PolePinchedSphere_p2.reg
)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
else(IFEM_COMMON_APP_BUILD)
  ifem_add_check_target()
endif()

# For generating the doxy
ifem_add_doc_target(TARGET ShellEl EXTRA_DIRS ${PROJECT_SOURCE_DIR}/..)

# Installation
include(GNUInstallDirs)
install(TARGETS ShellEl DESTINATION ${CMAKE_INSTALL_BINDIR})
