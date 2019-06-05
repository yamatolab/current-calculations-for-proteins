#! /bin/bash -e

rm -rf outdata
mkdir -p outdata

# calculate energy flux
time curp compute run.cfg > log

