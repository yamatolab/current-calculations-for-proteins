! /bin/bash

if [ "x$CURP_HOME" != "x" ]; then
    py_exec=$CURP_HOME/buildout/bin/curp-python
    ana_exec=$CURP_HOME/script/divide_flux.py
else
    py_exec=../../buildout/bin/curp-python
    ana_exec=../../script/divide_flux.py
fi

$py_exec  $ana_exec  -o outdata/flux_atm.dat -t 1.0 \
    outdata/flux_atm.dat00000  \
    > outdata/ana_atm.log

$py_exec  $ana_exec  -o outdata/flux_grp.dat -t 1.0 \
    outdata/flux_grp.dat00000  \
    > outdata/ana_grp.log
