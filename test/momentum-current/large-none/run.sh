#! /bin/bash -e

mkdir -p outdata
rm -f outdata/*.dat*

# calculate energy flux
time  curp compute  run.cfg > log

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

    curp analyze simplify-tensor \
        -i outdata/${grain}.dat > outdata/${grain}-sim.dat

done
