#! /bin/bash -e

script_dir="$CURP_HOME/script"
data_fp=outdata/grp-sim.dat

python $script_dir/plot_stress.py -i $data_fp -f outdata/grp-stress.pdf
python $script_dir/plot_ratio.py  -i $data_fp -f outdata/grp-ratio.pdf \
    -e 'total,kinetic'
