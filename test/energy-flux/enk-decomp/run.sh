#! /bin/bash -e

mkdir -p outdata
rm -f outdata/*.dat*

# calculate energy flux
time  curp compute  run.cfg > log

# postprocess
time  curp analyze divide-flux \
    -o outdata/flux.dat.gz -t 0.5 \
    outdata/flux_grp.dat00000  \
    > outdata/ana.log
