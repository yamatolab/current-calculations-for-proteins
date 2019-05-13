#! /bin/bash

$CURP_HOME/bin/ana-curp graph_ec.py -N -c \
    testdata/cluster.dat -g 1- -f testdata/ec.pdf \
    -n round ./testdata/easy-etc.dat; open ./testdata/ec.pdf
