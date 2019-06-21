#! /bin/bash

for dir in $(ls -d */)
do
    sed s/energy-flux/heat-flux/ ${dir}run.cfg > ${dir}hrun.cfg
    mv ${dir}run.cfg ${dir}erun.cfg
    sed -e 's/ run.cfg/ erun.cfg/' -i ${dir}run.sh
done
