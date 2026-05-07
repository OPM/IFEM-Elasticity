#!/bin/sh
# $Id$
##################################################################
# This script defines a set of sample simulation cases that
# is used to verify the integrity of an IFEM-based simulator.
# Copy this file to a sub-folder Test in your App-directory and
# insert the simulations that you want to use as regression tests.
##################################################################

# Define the name of the executable here
mysim=NonLinEl
if ! which $mysim; then mysim=../Release/bin/$mysim; fi

run () {
# This function runs a simulation with the specified options,
# pipes the terminal output to a log-file, and compares it
# with a previous simulation stored in the Reference folder.
  ext=`echo $1 | sed 's/.*\.//'`
  if `echo $2 | grep -q '^\-'`; then
    inp=$1
    log=`basename $1 .$ext`.log
    xinp=
  else
    inp=$2.$ext
    log=$2.log
    ln -s $1 $inp
    xinp=$inp
    shift
  fi
  echo Running $inp ...
  shift
  time -f "real time %E, user time %U, system time %S" \
  $mysim $inp $* > $log
  if [ $xinp ]; then rm $xinp; fi
  if [ ! -e Reference/$log ]; then
    mv $log Reference
  elif cmp -s $log Reference/$log; then
    echo Ok
  else
    echo "Warning: Discrepancies between current run and reference."
    diff $log Reference
  fi
}

if [ ! -d Reference ]; then mkdir Reference; fi

#####################################
# Enter the various cases below:
# Format: run <inputfile> [<options>]
#####################################

# 3D cantilever beam, linear-elastic material
for input in CanTS-p?.xinp; do
  ln -s $input TL-$input
  run TL-$input -dense -vtf 1
  ln -s $input UL-$input
  run UL-$input -UL -dense -vtf 1
  rm ?L-$input
done

# 3D Thick cylinder compression, hyperelastic material
run Cyl-p2.xinp -vtf 1 -nviz 3
run Cyl-p3.xinp -vtf 1 -nviz 3
run Cyl-p4.xinp -vtf 1 -nviz 3

# 2D Rubber block, hyperelastic material
# - mixed formulation with internal pressure modes
run FBlock-h9x5-p2.xinp FBlock-h9x5-Q2P1 -MX1
run FBlock-h8x3-p3.xinp FBlock-h8x3-Q3P2 -MX2
run FBlock-h8x2-p4.xinp FBlock-h8x2-Q4P3 -MX3
# - mixed formulation with continuous pressure field
run FBlock-h9x5-p2.xinp FBlock-h9x5-Q2Q1 -mixed
run FBlock-h8x3-p3.xinp FBlock-h8x3-Q3Q2 -mixed
run FBlock-h8x2-p4.xinp FBlock-h8x2-Q4Q3 -mixed

# Tension of a 2D elasto-plastic strip
run Necking-p2.xinp Necking-Q2P1 -MX1 -vtf 1 -nviz 3 -outPrec 6
run Necking-p2.xinp Necking-Q2Q1 -mixed -vtf 1 -nviz 3 -outPrec 6
run Necking-p2.xinp Necking-Q2-Q1 -Mixed -vtf 1 -nviz 3 -outPrec 6
run Necking-AxS-Fbar2.xinp -2Daxis -vtf 1 -nviz 3 -outPrec 6

# Cooks membrane, 2D plasticity
run Cook2D-p1-h1.xinp -vtf 1 -nviz 3 -outPrec 6
run Cook2D-p2-h1.xinp -vtf 1 -nviz 3 -outPrec 6
run Cook2D-h1-p2.xinp -vtf 1 -outPrec 6
