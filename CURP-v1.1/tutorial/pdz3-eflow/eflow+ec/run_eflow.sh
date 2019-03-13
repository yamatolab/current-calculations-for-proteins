#! /bin/bash -e

nproc=${1:-2}

mkdir -p outdata
rm -f outdata/*.dat*

# calculate energy flux
time mpiexec -n $nproc $CURP_HOME/bin/curp eflow.cfg > eflow.log

