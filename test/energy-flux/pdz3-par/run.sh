#! /bin/bash -e

nproc=${1:-2}

mkdir -p outdata
rm -f outdata/*.dat*

# calculate energy flux
# time curp compute run.cfg > ser.log
time mpiexec -n $nproc curp compute run.cfg > log

