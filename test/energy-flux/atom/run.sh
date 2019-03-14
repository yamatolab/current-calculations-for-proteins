#! /bin/bash -e

mkdir -p outdata
rm -f outdata/*.dat*

# calculate energy flux
time  $CURP_HOME/bin/curp  run.cfg > log

# postprocess
time  $CURP_HOME/bin/ana-curp  divide_flux.py \
    -o outdata/flux.dat.gz -t 1.0 \
    outdata/flux_atm.dat00000  \
    > outdata/ana.log
