#! /bin/bash

material_dir=$CURP_HOME/example/amber-b2AR

time  mpiexec -n 2 $CURP_HOME/bin/ana-curp  pickup_respairs.py \
    -p $material_dir/stripped.prmtop.gz -c 6.0 \
    $material_dir/nev.mdcrd.gz
