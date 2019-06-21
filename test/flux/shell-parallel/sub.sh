#! /bin/bash

# Usage:
#   $ sub.sh <nproc> [job system command and its options]

# Example:
#   $ sub.sh 8 bsub -o lava.out -m 'gtr10'

# if the submit of job failed, type the next command,
#   $ bkill $(seq 1559 1580)

crd_basefn=sam.mdcrd.gz
vel_basefn=sam.mdvel.gz

nproc=$1
shift
job_options=$@

for iproc in $(seq $nproc)
do
    crdtraj_fn=$iproc.$crd_basefn
    veltraj_fn=$iproc.$vel_basefn
    output=$iproc.outdata/flux.dat
    config_fn=$iproc.run.cfg

    cat run.cfg \
        | sed -e "s/{CRD_TRAJECTORY}/$crdtraj_fn/" \
        | sed -e "s/{VEL_TRAJECTORY}/$veltraj_fn/" \
        | sed -e "s%{OUTPUT}%$output%" \
        > $config_fn

    sleep 0.1

    $job_options  ./run.sh  $iproc
done

