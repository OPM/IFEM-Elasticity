# IFEM Elasticity


## Introduction

This module contains Elasticity libraries and applications built
using the IFEM library.

### Getting all dependencies

1. Install IFEM from https://github.com/OPM/IFEM

### Getting the code

This is done by first navigating to a folder `<App root>` in which you want
the application and typing

    git clone https://github.com/OPM/IFEM-Elasticity

### Compiling the code

To compile, first navigate to the root catalogue `<App root>`.

1. `cd IFEM-Elasticity/Linear`
2. `mkdir Debug`
3. `cd Debug`
5. `cmake -DCMAKE_BUILD_TYPE=Debug ..`
6. `make`

This will compile the libraries and the linear elasticity application.
The executable can be found in the 'bin' sub-folder.
Change all instances of `Debug` with `Release` to drop debug-symbols,
and get a faster running code.

### Testing the code

IFEM uses the cmake test system.
To compile and run all regression- and unit-tests, navigate to your build
folder (i.e. `<App root>/IFEM-Elasticity/Linear/Debug`) and type

    make check
