#! /usr/bin/env python
from __future__ import print_function

import os, sys
import numpy
import argparse

import curp_module

curp_srcdir = os.path.join(os.environ['CURP_HOME'], 'src')
if curp_srcdir not in sys.path:
    sys.path.insert(0, curp_srcdir)


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
from setproctitle import setproctitle
setproctitle('conv-trj')

def main():
    import os, sys

    # write the citation on log
    print()
    citation_fp = os.path.join(curp_srcdir,'..', 'LICENSE-short.txt')
    with open(citation_fp, 'rb') as citation_file:
        for line in citation_file:
            print('# '+line.strip())

    # parse arguments
    args = get_arguments()
    # print(args)

    # get topology and parameter file
    tpl = curp_module.get_topology(args.tpl_fmt, args.tpl_fp)

    # definition of the processing
    pname = args.proc_name
    if pname == 'dryrun':
        process = do_tryrun
    if pname == 'convert-only':
        process = do_convert_only
    if pname == 'mask':
        from mask import do_mask
        process = do_mask
    if pname == 'adjust-vel':
        from adjustvel import adjust_vel
        process = adjust_vel
    if pname == 'dist':
        from distance import do_distance
        process = do_distance

    # processing
    process(args, tpl)

def do_dryrun(args, tpl, trj=None):
    pass

def gen_trj(args, tpl):
    """Generate trajectory iterator."""

    if args.is_crd: trj_type = 'crd'
    if args.is_vel: trj_type = 'vel'

    trj_fns_list  = args.input_trj_fns
    fmts     = args.input_trj_fmts
    fst_lst_int_list = args.input_fst_lst_int

    # get number of atoms
    natom = tpl.get_natom()
    use_pbc = True if tpl.get_pbc_info() else False

    trj = curp_module.gen_trajectory_multi(
            trj_fns_list, fmts, natom, fst_lst_int_list, trj_type,
            use_pbc=use_pbc)

    return trj

def do_convert_only(args, tpl, trj=None):

    # get trajectory
    if trj is None:
        trj = gen_trj(args, tpl)

    # dt = trj.get_dt()
    dt = 0.01

    # write trajectory
    writer = curp_module.TrjWriter(args.output_trj_fn, args.output_trj_fmt, dt,
            args.is_vel, args.output_fst_lst_int)

    for ifrm, trj, box in trj:
        writer.write(ifrm-1, trj, box)
    writer.close()

def get_arguments():

    parser = argparse.ArgumentParser(description=(
        'Convert the trajectory into other trajectory that the CURP '+
        'can handle, with any processing.'))

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-crd',  # metavar='TRJ_TYPE',
            dest='is_crd', required=False, action='store_true',
            help=("specify the format of the trajectory."
                + "This argument allows you specify multiple formats."))

    group.add_argument('-vel',#  metavar='TRJ_TYPE',
            dest='is_vel', required=False, action='store_true',
            help=("specify the format of the trajectory."
                + "This argument allows you specify multiple formats."))

    # add argument definitions
    parser.add_argument('-i', '--input-filenames', metavar='TRJ_FILENAME',
            dest='input_trj_fns', required=True, nargs='+', action='append',
            help=('The trajectory file names.'))

    parser.add_argument('-if', '--input-formats', metavar='FORMAT',
            dest='input_trj_fmts', required=True, action='append',
            help=("specify the format of the trajectory."
                 +"This argument allows you specify multiple formats."))

    parser.add_argument('--irange', metavar=('FIRST','LAST','INTER'),
            dest='input_fst_lst_int', required=False,
            default=None, type=int, nargs=3,  action='append',
            help='The trajectory range to process over all trajectory file.')

    parser.add_argument('-o', '--output-filename', metavar='TRJ_FILENAME',
            dest='output_trj_fn', required=False,
            help=('The trajectory file name.'))

    parser.add_argument('-of', '--output-format', metavar='FORMAT',
            dest='output_trj_fmt', required=False,
            help=('The trajectory file format for output.'))

    parser.add_argument('--orange', metavar=('FIRST','LAST','INTER'),
            dest='output_fst_lst_int', required=False,
            default=[0,-1,1], type=int, nargs=3,
            help=('The trajectory range for output'))

    parser.add_argument('-pf', '--topology-format', metavar='FORMAT',
            dest='tpl_fmt', required=True,
            help='specify the format of the topology file.')

    parser.add_argument('-p', '--topology-file', metavar='TPL_FILENAME',
            dest='tpl_fp', required=True,
            help='specify the topology file.')
    
    # parser.add_argument('-m', '--method', metavar='METHOD',
            # dest='method', default='adjust-velocity',
            # required=False, choices=['adjust-velocity'],
            # help=('Processing method. The current version is available '
                # + 'for only adjusting velocity trajectory time.'))

    # adding of the processing command
    sp = parser.add_subparsers(dest='proc_name', help='processing command help')

    convert_only = sp.add_parser('convert-only',
            help='convert the trajectory into other format.')

    dryrun = sp.add_parser('dry-run',
            help='do dry-run mode.')

    adjust_vel = sp.add_parser('adjust-vel',
            help='Adjust time of the velocity trajectory to t from t-dt/2.')
    adjust_vel.add_argument('--any-arguments')

    mask = sp.add_parser('mask',
            help='Generate the trajectory that applies the mask.')
    mask.add_argument('-m', '--mask-filename', metavar='MASK_FILENAME',
            dest='mask_fn', required=True,
            help=('The mask file name.'))

    # distance
    dist = sp.add_parser('dist',
            help='Calculate inter-residue distances.')

    dist.add_argument('-m', '--method', metavar='DISTANCE_METHOD',
            dest='dist_method', required=True,
            choices=['cog', 'nearest', 'farthest'],
            help=('The method used in calculating the inter-resdue distances.'))

    dist.add_argument('-c', '--cutoff-length', metavar='CUTOFF_LENGTH',
            dest='dist_cutoff', required=True,
            default='5.0', type=float,
            help=('The distance cutoff in angstrom.'))

    parser.add_argument('-nf', '--name-format',
            dest='dist_format',
            required=False, default='{rid:05}_{rname}',
            help='specify the format for representing residue identify.')

    # make arguments
    return parser.parse_args()

def add_subcommand(sc_name, help, args):
    pass



if __name__ == '__main__':
    main()

    # test for adjust-vel subcommand

    # print('[if interval == 1]')
    # vels = ((i-1, 10*i, None) for i in range(1,20+1,1))
    # vel_pairs = gen_trj_late(vels, fst_lst_int=(5,-1,1))
    # for vel in vel_pairs:
        # print(vel)
    # print()

    # print('[if interval != 1]')
    # vels = ((i-1, 10*i, None) for i in range(1,20+1,1))
    # vel_pairs = gen_trj_late(vels, fst_lst_int=(2,-1,3))
    # for vel in vel_pairs:
        # print(vel)
    # print()


