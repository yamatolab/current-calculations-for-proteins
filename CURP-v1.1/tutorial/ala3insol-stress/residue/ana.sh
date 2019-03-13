#! /bin/bash -e
# postprocess

for filename in $(ls outdata/*.dat00000)
do
    dst_file=$(echo $filename | sed -e "s/00000//" | sed -e "s/stress_//")
    mv $filename $dst_file
done

for grain in atm inn out grp
do
    if [ ! -f outdata/${grain}.dat ]; then
        continue
    fi

    $CURP_HOME/bin/ana-curp simplify_tensor.py \
        -i outdata/${grain}.dat \
        -l 'bond,angle,torsion,improper,coulomb,vdw,kinetic' \
        outdata/${grain}_bond.dat     \
        outdata/${grain}_angle.dat    \
        outdata/${grain}_torsion.dat  \
        outdata/${grain}_improper.dat \
        outdata/${grain}_coulomb.dat  \
        outdata/${grain}_vdw.dat      \
        outdata/${grain}_kinetic.dat > outdata/${grain}-sim.dat

done
