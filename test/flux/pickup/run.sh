#! /bin/bash -e

material_dir=../../examples/amber-pdz3

time  curp analyze pickup-respairs \
    -p $material_dir/stripped.prmtop.gz -pf amber \
     -e 1:20 300:350 \
    -c 6.0 -if ascii $material_dir/strip.mdcrd.gz > 1.log

time  curp analyze pickup-respairs \
    -p $material_dir/stripped.prmtop.gz -pf amber \
     -e 1:20 300:350 -b \
    -c 6.0 -if ascii $material_dir/strip.mdcrd.gz > 2.log
