#! /bin/bash -e

example_path=../../examples/amber-enk-vacuum
prmtop_fp=$example_path/system.prmtop.gz

rm -rf outdata
mkdir -p outdata

# netcdf => netcdf
for trjtype in crd vel
do
    curp conv-trj -${trjtype} \
        -p $prmtop_fp  -pf amber \
        -i $example_path/sam.nc${trjtype} -if netcdf  --irange 0 -1 5 \
        -o outdata/reduce.${trjtype}.nc    -of netcdf  --orange 1 -1 1  \
        convert-only > ${trjtype}.ncdf.log

done

# netcdf => amber ascii
for trjtype in crd vel
do
    curp conv-trj -${trjtype} \
        -p $prmtop_fp  -pf amber \
        -i $example_path/sam.nc${trjtype} -if netcdf  --irange 0 -1 5 \
        -o outdata/reduce.md${trjtype}    -of ascii   --orange 1 -1 1  \
        convert-only > ${trjtype}.ascii.log

done

