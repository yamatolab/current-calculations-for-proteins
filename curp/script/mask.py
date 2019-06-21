#! /usr/bin/env python
from __future__ import print_function

import os, sys
from curp import TrjWriter


def do_mask(tpl, trj, trj_type, output_trj_fn,
            output_trj_fmt, mask_fn, output_fst_lst_int=[(0,-1,1)], **kwds):

    # dt = trj.get_dt()
    dt = 0.01
    is_vel = trj_type=='vel'
    ids = load_mask(mask_fn)

    # print(ids)

    # apply the mask to the trajectory
    # def gen_mask_trj_and_print(trj):
        # for istp, snap, box in trj:
            # print('istp = ', istp)
    trj_new = ( (istp,snap[ids],box) for (istp,snap,box) in trj )

    # write trajectory
    writer = TrjWriter(output_trj_fn, output_trj_fmt, dt,
            is_vel, output_fst_lst_int)

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
