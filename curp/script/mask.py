#! /usr/bin/env python
from __future__ import print_function

import os, sys
import curp_module

curp_srcdir = os.path.join(os.environ['CURP_HOME'], 'src')
if curp_srcdir not in sys.path:
    sys.path.insert(0, curp_srcdir)


def do_mask(args, tpl, trj=None):

    if trj is None:
        import conv_trj
        trj = conv_trj.gen_trj(args, tpl)

    # dt = trj.get_dt()
    dt = 0.01

    fn = args.mask_fn 
    ids = load_mask(fn)

    # print(ids)

    # apply the mask to the trajectory
    # def gen_mask_trj_and_print(trj):
        # for istp, snap, box in trj:
            # print('istp = ', istp)
    trj_new = ( (istp,snap[ids],box) for (istp,snap,box) in trj )

    # write trajectory
    writer = curp_module.TrjWriter(args.output_trj_fn, args.output_trj_fmt, dt,
            args.is_vel, args.output_fst_lst_int)

    for ifrm, trj, box in trj_new:
        writer.write(ifrm-1, trj, box)
    writer.close()

def load_mask(fn):

    with open(fn, 'rb') as file:
        ids_str = (line.split()[1] for line in file
                if line.startswith('ATOM'))

        ids = [ int(id_str)-1 for id_str in ids_str ]

    return ids

if __name__ == '__main__':
    pass
