#! /bin/bash -e

mkdir -p outdata
rm -f outdata/*.dat*

# output config file with default values
$CURP_HOME/bin/curp  --output-conf-default > outdata/default.cfg

# output config file with default values in reStructuredText format
$CURP_HOME/bin/curp  --output-conf-formatted > outdata/parameters.rst
