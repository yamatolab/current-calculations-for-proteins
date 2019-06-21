#! /bin/bash

trap 'exit 1' INT


if [ -z "$CURP_HOME" ]; then
    echo "Please set CURP_HOME environment variable."
    exit 1
fi

for category in momentum-current  flux  dynamics  util
do
    cd $category
    echo "** entering $category **"
    ./runall.sh
    cd ../
done
