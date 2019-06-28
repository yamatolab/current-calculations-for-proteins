#! /bin/bash

trap 'exit 1' INT

for category in momentum-current  flux  dynamics  util
do
    cd $category
    echo "** entering $category **"
    ./runall.sh
    cd ../
done
