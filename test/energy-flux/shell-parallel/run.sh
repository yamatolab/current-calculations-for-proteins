#! /bin/bash -e

iproc=$1

echo $iproc

rm -rf $iproc.outdata

if [ ! -d $iproc.outdata ]; then
    mkdir $iproc.outdata
fi
time $CURP_HOME/bin/curp  -s  $iproc.run.cfg > $iproc.log
