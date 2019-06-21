#! /bin/bash -e

mkdir -p outdata
rm -f outdata/*.dat*

# output config file with default values
curp compute  --output-conf-default > outdata/default.cfg

# output config file with default values in reStructuredText format
curp compute  --output-conf-formatted > outdata/parameters.rst
