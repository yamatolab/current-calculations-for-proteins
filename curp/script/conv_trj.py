#! /usr/bin/env python
from __future__ import print_function

import os, sys
import numpy
import argparse

from curp import get_tpl, TrjWriter
import curp.parser as parser


"""
Supported parameter and topology format:
    - pdb
    - amber

Supported trajectory format:
    - restart: amber restart
    - ascii: amber trajectory in ascii format
    - netcdf: amber trajectory in netcdf format
"""


# make the process name the CURP.

def do_dryrun(tpl, trj=None, **kwds):
    pass

def gen_trj(tpl, input_trj_fns, input_trj_fmts, trj_type,
            input_fst_lst_int=[(1, -1, 1)], use_pbc=True, **kwds):
    """Generate trajectory iterator"""

    trj_fns_list  = input_trj_fns
    fmts     = input_trj_fmts
    fst_lst_int_list = input_fst_lst_int

    # get number of atoms
    natom = tpl.get_natom()
    use_pbc = True if (tpl.get_pbc_info() and use_pbc) else False

    trj = gen_trajectory_multi(
            trj_fns_list, fmts, natom, fst_lst_int_list, trj_type,
            use_pbc=use_pbc)

    return trj

def gen_trajectory_multi(trj_fns_list, trj_fmts, natom, fst_lst_int_list,
        trj_type='crd',logger=None, use_pbc=False):

    from itertools import izip_longest

    # print(" trj_fns: ", trj_fns_list,
    #       "\n trj_fmts: ", trj_fmts,
    #       "\n fst_lst_int_list: ", fst_lst_int_list)
    for trj_fns, fmt, fst_lst_int in izip_longest(
            trj_fns_list, trj_fmts, fst_lst_int_list):

        # determine parser
        Parser = get_any_parser(fmt, natom, use_pbc, trj_type)

        # generate trajectory
        trj = parser.gen_trajectory(trj_fns, Parser, natom, logger=None,
                use_pbc=use_pbc, first_last_interval=fst_lst_int)

        # iterate
        for istp, snap, box in trj:
            yield istp, snap, box

def get_any_parser(trj_fmt, natom, use_pbc, trj_type='crd'):

    if trj_fmt == 'restart':
        Parser = parser.get_parser(trj_fmt, 'restart')
        Parser.set_trjtype(trj_type)
        return Parser

    else:
        if trj_type == 'crd':
            return parser.get_parser(trj_fmt, 'coordinate')

        elif trj_type == 'vel':
            return parser.get_parser(trj_fmt, 'velocity')

def do_convert_only(tpl, trj, output_trj_fn, output_trj_fmt,
                    trj_type, output_fst_lst_int=[0,-1,1], **kwds):

    # dt = trj.get_dt()
    dt = 0.01
    is_vel = trj_type == 'vel'

    # write trajectory
    writer = TrjWriter(output_trj_fn, output_trj_fmt, dt,
                                   is_vel, output_fst_lst_int)

    for ifrm, trj, box in trj:
        writer.write(ifrm-1, trj, box)
    writer.close()


if __name__ == '__main__':
    from console import arg_conv_trj, exec_command

    parser = arg_conv_trj
    exec_command(parser)
