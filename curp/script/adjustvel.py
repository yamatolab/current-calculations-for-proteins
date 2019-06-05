#! /usr/bin/env python
from __future__ import print_function

import os, sys
from curp import TrjWriter

def adjust_vel(tpl, trj, trj_type, output_trj_fn,
               output_trj_fmt, output_fst_lst_int=[(0,-1,1)], **kwds):

    # get velocity trajectory
    is_vel = trj_type=='vel'
    if not is_vel:
        msg = "Adjust velocity command doesn't support coordinate trajectory."
        raise Exception(msg)

    # dt = trj.get_dt()
    dt = 0.01

    # process velocity trajectory
    vel_pair_iter = gen_trj_late(trj, output_fst_lst_int)

    def gen_adjust_vel(vel_pair_iter):
        for trj_prev, trj_now in vel_pair_iter:
            istp_p, vel_p, box_p = trj_prev
            istp,   vel,   box   = trj_now
            # if box_p:
                # box = 0.5*(box+box_p)

            print('adjust-vel: ', istp_p, istp)
            yield istp_p, (0.5*(vel_p + vel), box)

    adjusted_vels = gen_adjust_vel(vel_pair_iter)

    # write trajectory
    writer = TrjWriter(output_trj_fn, output_trj_fmt, dt,
                                   is_vel, (1,-1,1))

    for ifrm, (vel, box) in adjusted_vels:
        writer.write(ifrm-1, vel)
    writer.close()

def gen_trj_late(trj, fst_lst_int=(1,-1,1)):

    first, last, inter = fst_lst_int

    nstep = 0
    first = 0 if first < 0 else first
    last  = 10**18 if last <= -1 else last

    if inter==1:
        for istp_cur, snap, box in trj:
            nstep += 1

            if nstep < first: continue

            if nstep > last:
                yield (istp_prev, snap_prev, box_prev), (istp_cur, snap, box)
                break

            if nstep == first:
                istp_prev = istp_cur
                snap_prev = snap
                box_prev = box
                continue

            yield (istp_prev, snap_prev, box_prev), (istp_cur, snap, box)
            istp_prev = istp_cur
            snap_prev = snap
            box_prev = box

    else:
        cnt_inter = 0
        for istp_cur, snap, box in trj:
            nstep += 1

            if nstep < first: continue
            if nstep > last: break

            cnt_inter += 1
            if nstep == first or cnt_inter == inter:
                cnt_inter = 0
                istp_prev = istp_cur
                snap_prev = snap
                box_prev = box

                # one more step
                nstep += 1
                cnt_inter += 1
                istp_cur, snap, box = trj.next()
                yield (istp_prev, snap_prev, box_prev), (istp_cur, snap, box)

if __name__ == '__main__':

    # test for adjust-vel subcommand
    print('[if interval == 1]')
    vels = ((i-1, 10*i, None) for i in range(1,20+1,1))
    vel_pairs = gen_trj_late(vels, fst_lst_int=(5,-1,1))
    for vel in vel_pairs:
        print(vel)
    print()

    print('[if interval != 1]')
    vels = ((i-1, 10*i, None) for i in range(1,20+1,1))
    vel_pairs = gen_trj_late(vels, fst_lst_int=(2,-1,3))
    for vel in vel_pairs:
        print(vel)
    print()


