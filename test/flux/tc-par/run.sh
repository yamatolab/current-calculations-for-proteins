#! /bin/bash -e

nproc=${1:-2}

mkdir -p outdata
rm -f outdata/*.dat*

time mpirun -n $nproc  curp cal-tc --no-axes\
    --frame-range 1 50000 1 \
    --average-shift 1 --sample-number 0 \
    -o outdata/tc.dat -a outdata/acf.nc flux_grp.nc > log
