#! /bin/bash -e

nproc=${1:-2}

# postprocess
time mpiexec -n $nproc $CURP_HOME/bin/cal-tc \
    --frame-range 1 50 1 --average-shift 1 \
    -a outdata/acf.nc \
    -o outdata/ec.dat outdata/flux_grp.nc > ec.log
