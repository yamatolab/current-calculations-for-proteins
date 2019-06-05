#! /bin/bash -e

mkdir -p outdata
rm -f outdata/*.dat*

# calculate energy flux
time  curp compute  run.cfg > log

# postprocess
for grain in atm grp
do
    time  curp analyze divide-flux \
        -o outdata/flux_${grain}.dat.gz -t 1.0 \
        outdata/flux_${grain}.dat00000  \
        > outdata/ana_${grain}.log
done
