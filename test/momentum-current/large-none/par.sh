#! /bin/bash

mkdir -p outdata
rm -f outdata/*.dat*

# calculate energy flux
time  mpiexec -n 2 $CURP_HOME/bin/curp  run.cfg > log

# postprocess
for filename in $(ls outdata/*.dat00000)
do
    dst_file=$(echo $filename | sed -e "s/00000//" | sed -e "s/current_//")
    mv $filename $dst_file
done

for grain in atm grp
do
    if [ ! -f outdata/${grain}.dat ]; then
        continue
    fi

    $CURP_HOME/bin/ana-curp simplify_tensor.py \
        -i outdata/${grain}.dat > outdata/${grain}-sim.dat

done
