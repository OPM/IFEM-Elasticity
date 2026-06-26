## Introduction

This repository contains applications and libraries for solving various solid
mechanics problems, built using the [IFEM](https://github.com/OPM/IFEM) library.

### Getting all dependencies

1. Install IFEM from https://github.com/OPM/IFEM

### Getting the code

This is done by first navigating to a folder `<App root>` in which you want
the applications and typing

    git clone https://github.com/OPM/IFEM-Elasticity

### Compiling the code

To compile, first navigate to the root catalogue `<App root>`.

    cd IFEM-Elasticity
    mkdir Debug
    cd Debug
    cmake -DCMAKE_BUILD_TYPE=Debug ..
    make

This will compile the libraries and three applications.
The executables can be found in the `bin` sub-folder.
Change all instances of `Debug` with `Release` to drop debug-symbols,
and get a faster running code.

### Testing the code

IFEM uses the cmake test system.
To compile and run all regression- and unit-tests, navigate to your build
folder (i.e., `<App root>/IFEM-Elasticity/Debug`) and type

    make check
