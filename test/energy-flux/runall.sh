#! /bin/bash

red="\033[31m" # red
gre="\033[32m" # green
off="\033[0m"  # color off

trap 'rm -rf /tmp/curptest.*; exit 1' INT

for dir in $(echo *)
do
    [ "$dir" == "shell-parallel" ] && continue
    if [ -d $dir ]; then
        echo -n "performing $dir ..."
        cd $dir
        error_fp=$(mktemp /tmp/curptest.${dir}.XXXXXX)

        # run
        ./run.sh > /dev/null 2> $error_fp

        # error check
        error=$?
        if [ $error -eq 0 ]; then
            echo -e "\r\c"
            echo -e "[${gre}PASSED${off}] $dir                  "
        else
            echo ""
            cat $error_fp
            echo -e "[${red}FAILED${off}] $dir                  "
        fi
        rm -f $error_fp

        cd ../
    fi
done
