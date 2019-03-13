#! /bin/bash -e

mkdir -p outdata
rm -f outdata/*.dat*

# calculate energy flux
time  $CURP_HOME/bin/curp  run.cfg > log
