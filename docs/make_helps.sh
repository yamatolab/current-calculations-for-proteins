#! /bin/bash

help_text=source/helps.txt

rm -f $help_text
for cmd in "$@"
do
    echo $cmd
    cmd_line=$(basename $cmd)
    cmd_line=${cmd_line/.py/}
    echo $cmd_line         >> $help_text
    echo '---------------' >> $help_text
    echo ''       >> $help_text
    echo ''       >> $help_text
    echo '::'     >> $help_text
    echo ''       >> $help_text

    if [ "${cmd##*.}" == "py" ]; then
        $CURP_HOME/bin/ana-curp ${cmd} --help \
            | awk '{printf"   %s\n",$0}' >> $help_text
    else
        ${cmd} --help | awk '{printf"   %s\n",$0}' >> $help_text
    fi

    echo ''       >> $help_text
    echo ''       >> $help_text
done
