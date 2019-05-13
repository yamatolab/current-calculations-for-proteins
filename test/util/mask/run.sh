#! /bin/bash -e

example_path=$CURP_HOME/test/amber-enk-vacuum
prmtop_fp=$example_path/system.prmtop.gz

rm -rf outdata
mkdir -p outdata

$CURP_HOME/bin/conv-trj -crd \
    -p $prmtop_fp  -pf amber \
    -i $example_path/sam.nccrd -if netcdf  --irange 1 -1 1 \
    -o outdata/masked.crd.nc   -of netcdf  --orange 1 -1 1 \
    mask -m mask.pdb > 1.log

$CURP_HOME/bin/conv-trj -vel \
    -p $prmtop_fp  -pf amber \
    -i $example_path/sam.ncvel -if netcdf  --irange 1 -1 1 \
    -o outdata/masked1.vel.nc  -of netcdf  --orange 0 -1 5  \
    mask -m mask.pdb > 2.log
    
# with multiple format for input
$CURP_HOME/bin/conv-trj -vel \
    -p $prmtop_fp  -pf amber \
    -i $example_path/pre.rst   -if restart  --irange 1 -1 1 \
    -i $example_path/sam.ncvel -if netcdf   --irange 1 -1 1 \
    -o outdata/masked2.vel.nc  -of netcdf  --orange 1 -1 5  \
    mask -m mask.pdb > multi.log
