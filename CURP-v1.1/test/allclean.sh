#! /bin/bash

for category in energy-flux momentum-current util dynamics
do
    cd $category
    ./allclean.sh
    cd ..
done
