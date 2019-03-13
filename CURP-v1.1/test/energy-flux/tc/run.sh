#! /bin/bash -e

mkdir -p outdata
rm -f outdata/*.dat*

# calculate energy flux
time  $CURP_HOME/bin/curp  run.cfg > log

# postprocess
time  $CURP_HOME/bin/cal-tc \
    --frame-range 1 10 1 --average-shift 1 --sample-number 0 \
    -a outdata/acf0.nc \
    -o outdata/tc0.dat outdata/flux_grp.nc > tc0.log

# # postprocess
time  $CURP_HOME/bin/cal-tc \
    --frame-range 1 10 1 --average-shift 2 --sample-number 0 \
    -a outdata/acf1.nc \
    -o outdata/tc1.dat outdata/flux_grp.nc > tc1.log

time  $CURP_HOME/bin/ana-curp  sum_tc.py \
    outdata/tc0.dat outdata/tc1.dat > outdata/tc.dat

time $CURP_HOME/bin/sum-acf \
    -a outdata/acf.nc -t outdata/tcs.nc \
    outdata/acf0.nc outdata/acf1.nc

time $CURP_HOME/bin/ana-curp get_ncdata.py \
    outdata/acf.nc -n acf -r 1:2 \
    -o outdata/acf

time $CURP_HOME/bin/ana-curp get_ncdata.py \
    outdata/tcs.nc -n tcs -r 2:3 \
    -o outdata/tcs
