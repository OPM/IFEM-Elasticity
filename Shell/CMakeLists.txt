cmake_minimum_required(VERSION 3.22)

project(ShellEl LANGUAGES C CXX)

ifem_add_library(
  NAME
    Shell
  SOURCES
    KirchhoffLoveShell.C
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
