#! /bin/bash -e


for dir in chart-een-apo chart-een-a3rem chart-een-diff view-een-3D
do
    cd $dir
    make clean
    cd ../
done
