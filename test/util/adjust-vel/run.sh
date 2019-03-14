#! /bin/bash -e

example_path=$CURP_HOME/test/amber-enk-vacuum
prmtop_fp=$example_path/system.prmtop.gz
vel_fns=$example_path/sam.ncvel
output_fn=outdata/adjusted.vel.nc

mkdir -p outdata

$CURP_HOME/bin/conv-trj -vel \
    -p $prmtop_fp  -pf amber \
    -i $example_path/pre.rst   -if restart  --irange 1 -1 1 \
    -i $example_path/sam.ncvel -if netcdf   --irange 1 -1 1 \
    -o $output_fn  -of netcdf  --orange 1 -1 5  \
    adjust-vel > log
    
