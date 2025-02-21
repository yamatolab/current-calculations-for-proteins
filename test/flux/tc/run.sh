#! /bin/bash -e

mkdir -p outdata
rm -f outdata/*.dat*

# calculate energy flux
time  curp compute  erun.cfg > elog

# calculate heat flux
time  curp compute  hrun.cfg > hlog

# postprocess
time  curp cal-tc \
    --frame-range 1 10 1 --average-shift 1 --sample-number 0 \
    -a outdata/acf0.nc \
    -o outdata/tc0.dat outdata/eflux_grp.nc > tc0.log

# # postprocess
time  curp cal-tc \
    --frame-range 1 10 1 --average-shift 2 --sample-number 0 \
    -a outdata/acf1.nc \
    -o outdata/tc1.dat outdata/eflux_grp.nc > tc1.log

time curp analyze sum-tc outdata/tc0.dat outdata/tc1.dat > outdata/tc.dat

time curp sum-acf \
    -a outdata/acf.nc -t outdata/tcs.nc \
    outdata/acf0.nc outdata/acf1.nc

time curp analyze get-ncdata \
    outdata/acf.nc -n acf -r 1:2 \
    -o outdata/acf

time curp analyze get-ncdata \
    outdata/tcs.nc -n tcs -r 2:3 \
    -o outdata/tcs
