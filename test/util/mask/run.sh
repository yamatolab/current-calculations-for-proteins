#! /bin/bash -e

example_path=../../examples/amber-enk-vacuum
prmtop_fp=$example_path/system.prmtop.gz

rm -rf outdata
mkdir -p outdata

curp conv-trj -crd \
    -p $prmtop_fp  -pf amber \
    -i $example_path/sam.nccrd -if netcdf  --irange 1 -1 1 \
    -o outdata/masked.crd.nc   -of netcdf  --orange 1 -1 1 \
    mask -m mask.pdb > 1.log

curp conv-trj -vel \
    -p $prmtop_fp  -pf amber \
    -i $example_path/sam.ncvel -if netcdf  --irange 1 -1 1 \
    -o outdata/masked1.vel.nc  -of netcdf  --orange 0 -1 5  \
    mask -m mask.pdb > 2.log
    
# with multiple format for input
curp conv-trj -vel \
    -p $prmtop_fp  -pf amber \
    -i $example_path/pre.rst   -if restart  --irange 1 -1 1 \
    -i $example_path/sam.ncvel -if netcdf   --irange 1 -1 1 \
    -o outdata/masked2.vel.nc  -of netcdf  --orange 1 -1 5  \
    mask -m mask.pdb > multi.log
