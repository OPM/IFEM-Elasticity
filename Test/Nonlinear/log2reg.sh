#!/bin/sh
# $Id$
# This script creates new reg-files from specified log-files,
# for the FiniteDeformation simulations.

for file in $*; do if [ -f $file ]; then found=yes
echo $file
sed -f ../../../../../scripts/log2reg.sed $file | sed \
'1 s/^.*NonLinEl //;1 s/ *-vtf 1//;1 s/ *-nviz [1-9]//;1 a\
' > `basename $file .log`.reg
fi; done
if [ -z "$found" ]; then echo "usage $0 <logfile(s)>"; fi
