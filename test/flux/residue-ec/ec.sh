#! /bin/bash

for eflux_fn in outdata/*.flux.dat.gz
do
    
    don_acc_line=$(basename $eflux_fn)
    don_acc_line=${don_acc_line/.flux.dat.gz}

    $CURP_HOME/bin/ana-curp  cal_tc.py  $tc_exec \
        --dt 0.01  --frame-range 1 10 1 \
        --shift 1 --sample-number 0 \
        -a outdata/acf_${don_acc_line}.dat \
        -t outdata/tcs_${don_acc_line}.dat \
        -o outdata/log_${don_acc_line}.log   \
        $eflux_fn > outdata/tc_${don_acc_line}.dat
done

cat outdata/tc_*.dat | grep -v '%' > outdata/tc_all.dat

