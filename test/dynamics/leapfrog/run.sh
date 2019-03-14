#! /bin/bash -e

rm -rf outdata
mkdir -p outdata

# calculate energy flux
time  $CURP_HOME/bin/curp  run.cfg > log

