#! /bin/bash

mkdir -p outdata
rm -f outdata/*.dat*

# calculate energy flux
time  mpiexec -n 2 $CURP_HOME/bin/curp  run.cfg > log

# postprocess
for grain in atm grp
do
    time  $CURP_HOME/bin/ana-curp  divide_flux.py \
        -o outdata/flux_${grain}.dat.gz -t 1.0 \
        outdata/flux_${grain}.dat00000  \
        > outdata/ana_${grain}.log
done
