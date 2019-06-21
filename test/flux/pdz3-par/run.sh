#! /bin/bash -e

nproc=${1:-2}

mkdir -p outdata
rm -f outdata/*.dat*

# calculate energy flux
# time curp compute erun.cfg > ser.eelog

# calculate heat flux
# time curp compute hrun.cfg > ser.eelog
time mpiexec -n $nproc curp compute erun.cfg > elog

# calculate heat flux
# time curp compute hrun.cfg > ser.elog
time mpiexec -n $nproc curp compute hrun.cfg > hlog


