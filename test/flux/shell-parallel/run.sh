#! /bin/bash -e

iproc=$1

echo $iproc

rm -rf $iproc.outdata

if [ ! -d $iproc.outdata ]; then
    mkdir $iproc.outdata
fi
time curp compute  -s  $iproc.run.cfg > $iproc.log
