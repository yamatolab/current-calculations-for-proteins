#! /bin/bash

for dir in $(echo *)
do
    if [ -d $dir ]; then
        cd $dir
        echo $dir
        rm -rf log outdata
        cd ../
    fi
done
