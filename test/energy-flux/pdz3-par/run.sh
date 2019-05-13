#! /bin/bash -e

nproc=${1:-2}

mkdir -p outdata
rm -f outdata/*.dat*

# calculate energy flux
# time $CURP_HOME/bin/curp run.cfg > ser.log
time mpiexec -n $nproc $CURP_HOME/bin/curp run.cfg > log

