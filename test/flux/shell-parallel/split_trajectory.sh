#! /bin/bash

# Usage:
#    $ ./split_trajectory.sh <system.prmtop>  <mdcrd or mdvel> <first and last step of trajectory> <number of process for parallel calculation>

# Example:
#    $ ./split_trajectory.sh ../../amber-large/system.prmtop ../../amber-large/sam.mdcrd.gz 1,1000 15 > log

#    $ ./split_trajectory.sh ../../amber-large/system.prmtop ../../amber-large/sam.mdvel.gz 1,1000 15 > log


tpl_file=$1   # ../../system.prmtop
traj_file=$2  # ../../sam.mdcrd.gz
traj_range=$3 # 1,1000'
nproc=$4      # 50

box_arg=''
if echo $traj_file | grep 'vel'
then
    box_arg=nobox
fi

step_fst=$(echo $traj_range | cut -d , -f 1)
step_lst=$(echo $traj_range | cut -d , -f 2)

# make all range list
ntraj=$(($step_lst - $step_fst + 1))
ntraj_per_proc_1=$(($ntraj / $nproc))
ntraj_last_proc=$(($ntraj % $nproc))

#
traj_basefn=$(basename $traj_file)

for iproc in $(seq $nproc)
do
    step_beg=$(($ntraj_per_proc_1 * ($iproc - 1) + 1))
    step_end=$(($ntraj_per_proc_1 * $iproc))
    trajout_fn=$iproc.$traj_basefn

$AMBERHOME/bin/ptraj  $tpl_file  <<END_PTRAJ
trajin $traj_file  $step_beg  $step_end
trajout $trajout_fn  $box_arg
go
END_PTRAJ

done

if [ $ntraj_last_proc != '0' ]; then

    step_beg=$(($ntraj_per_proc_1 * $iproc + 1))
    step_end=$step_lst
    trajout_fn=$(($iproc+1)).$traj_basefn

$AMBERHOME/bin/ptraj  $tpl_file  <<END_PTRAJ
trajin $traj_file  $step_beg  $step_end
trajout $trajout_fn  $box_arg
go
END_PTRAJ

fi

