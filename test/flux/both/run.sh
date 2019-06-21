#! /bin/bash -e

mkdir -p outdata
rm -f outdata/*.dat*

# calculate energy flux
time  curp compute  erun.cfg > elog

# calculate heat flux
time  curp compute  hrun.cfg > hlog

# postprocess
for grain in atm grp
do
    time  curp analyze divide-flux \
        -o outdata/eflux_${grain}.dat.gz -t 1.0 \
        outdata/eflux_${grain}.dat00000  \
        > outdata/ana_${grain}.log
done
