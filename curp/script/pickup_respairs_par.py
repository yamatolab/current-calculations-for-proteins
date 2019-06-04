#! /usr/bin/env python
from __future__ import print_function

import os, sys
import numpy
import argparse

import curp

curp_srcdir = os.path.join(os.environ['CURP_HOME'], 'src')
if curp_srcdir not in sys.path:
    sys.path.insert(0, curp_srcdir)
from table.group import gen_residue_group

def main():
    import os, sys

    # parse arguments
    args = get_arguments()

    # get parameters
    cutoff = float(args.cutoff)
    cutoff2 = cutoff*cutoff

    # get topology
    tpl = curp.get_tpl(args.prmtop_fn)
    natom = tpl.get_natom()

    # generate residue group
    res_info  = tpl.get_residue_info()
    atom_info = tpl.get_atom_info()
    rname_iatoms_pairs = list( (rname, list(numpy.array(iatoms)-1))
            for rname, iatoms in gen_residue_group(
                res_info, atom_info, name_fmt=args.format))

    resnames = [ resname for resname, iatoms in rname_iatoms_pairs ]

    # get trajectory
    crd_parser = curp.gen_trj(
            args.traj_fns, natom, args.interval, use_pbc=False)
    crds = ( crd for itraj, crd, box in crd_parser )

    # generate residue pair table
    import lib_pickup
    cutoff_method = dict(
            com      = lib_pickup.is_com,
            nearest  = lib_pickup.is_nearest,
            farthest = lib_pickup.is_farthest)[args.cutoff_method]

    import parallel
    par = parallel.ParallelProcessor()

    method = gen_respair_table_par(
            resnames, rname_iatoms_pairs, cutoff_method, cutoff2)
    respair_table_traj = par.run(method, crd=crds)

    if par.is_root():
        respair_iter = gen_respair_table_over_traj(
                respair_table_traj,args.is_union)

        # write residue pair data
        print('# residue pair within {:4.2f}'.format(cutoff))
        npair = 0
        for rname_i, rnames_j in respair_iter:
            if rnames_j:
                npair += len(rnames_j)
                print('[{}]'.format(rname_i))
                print('  '.join(rnames_j))

        print('# number of pairs = ', npair)

    single = False
    if single:

        respair_table_traj = ( gen_respair_table( crd, resnames,
                    rname_iatoms_pairs, cutoff_method, cutoff2)
                    for itraj, crd, box in crd_parser )

        respair_iter = gen_respair_table_over_traj(
                respair_table_traj, args.is_union)


        # write residue pair data
        print('# residue pair within {:4.2f}'.format(cutoff))
        npair = 0
        for rname_i, rnames_j in respair_iter:
            if rnames_j:
                npair += len(rnames_j)
                print('[{}]'.format(rname_i))
                print('  '.join(rnames_j))

        print('# number of pairs = ', npair)


def get_arguments():

    parser = argparse.ArgumentParser(
            ('Pick up the residue pair table to given file as a ndx format.'))

    # add argument definitions
    parser.add_argument(
            dest='traj_fns', nargs='+',
            help='specify the trajectory files.')

    parser.add_argument(
            '-p', '--input-prmtop-file', dest='prmtop_fn', required=True,
            help='specify the prmtop file to be amber format.')
    
    parser.add_argument(
            '-i', '--interval', dest='interval', default=1, type=int,
            help=('specify interval step to perform the calculation for '
                  'trajectory.'))

    parser.add_argument(
            '-m', '--cutoff-method', dest='cutoff_method', default='nearest',
            required=False, choices=['com', 'nearest', 'farthest'],
            help='cutoff method; com, nearest and farthest for residue.')

    parser.add_argument(
            '-c', '--cutoff', dest='cutoff',
            required=False, default=5.0,
            help='specify the cutoff to pick up.')

    parser.add_argument(
            '-t', '--trim-resnames', dest='trim_resnames', nargs='*', 
            required=False, default=[],
            help='residue names for trimming from the target residues.')

    parser.add_argument(
            '-f', '--format', dest='format',
            required=False, default='{rid:05}_{rname}',
            help='specify the format for representing residue identify.')

    parser.add_argument(
            '-U', '--disble-union', dest='is_union',
            required=False, action='store_false',
            help=('whether residue pair table is calculated by union '
                  'or intersection when collecting over trajectoryspecify'))

    # make arguments
    return parser.parse_args()

def cal_dist2(p1, p2):
    # return sum( (p1 - p2)**2 )
    return (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2

def zip_double(data_iter):

    for i, e1 in enumerate(data_iter):
        for j, e2 in enumerate(data_iter):
            if j <= i: continue
            yield e1, e2

def gen_respair_table(crd, resnames, rname_iatoms_pairs, cutoff_method,cutoff2):

    res_crds = [ crd[iatoms] for rname, iatoms in rname_iatoms_pairs ]

    for i, (rname_i, rcrd_i) in enumerate(zip(resnames, res_crds)):
        resj_cands = []
        for j, (rname_j, rcrd_j) in enumerate(zip(resnames, res_crds)):
            if j <= i: continue

            if cutoff_method(rcrd_i, rcrd_j, cutoff2):
                resj_cands.append(rname_j)

        yield rname_i, resj_cands

def gen_respair_table_par(resnames, rname_iatoms_pairs, cutoff_method, cutoff2):

    def wrapper(crd):
        res_crds = [ crd[iatoms] for rname, iatoms in rname_iatoms_pairs ]
        table = []

        for i, (rname_i, rcrd_i) in enumerate(zip(resnames, res_crds)):
            resj_cands = []
            for j, (rname_j, rcrd_j) in enumerate(zip(resnames, res_crds)):
                if j <= i: continue

                if cutoff_method(rcrd_i, rcrd_j, cutoff2):
                    resj_cands.append(rname_j)

            table += [(rname_i, resj_cands)]

    return wrapper

def is_com(cutoff2):

    def _judge(crd_i, crd_j):
        com_i = numpy.sum(crd_i, 0)
        com_j = numpy.sum(crd_j, 0)
        return cal_dist2(com_i, com_j) <= cutoff2

    return _judge

def is_nearest(cutoff2):

    def _judge(crd_i, crd_j):
        for r_i in crd_i:
            for r_j in crd_j:
                if cal_dist2(r_i, r_j) <= cutoff2:
                    return True
            
        else:
            return False

    return _judge

def is_farthest(cutoff2):

    def _judge(crd_i, crd_j):
        for r_i in crd_i:
            for r_j in crd_j:
                if cal_dist2(r_i, r_j) >= cutoff2:
                    return False
                
        else:
            return True

    return _judge

def gen_respair_table_over_traj(respair_table_traj, is_union=True):

    table = {}
    rnames = []

    for itrj_1, respair_table in enumerate(respair_table_traj):
        print('# itraj = ', itrj_1+1)

        for rname_i, rnames_j in respair_table:

            if rname_i not in table:
                table[rname_i] = set(rnames_j)
                rnames.append( rname_i )
                continue

            if is_union:
                table[rname_i] = table[rname_i] | set(rnames_j)
            else:
                table[rname_i] = table[rname_i] & set(rnames_j)

    for rname_i in rnames:
        yield rname_i, list(table[rname_i])

if __name__ == '__main__':
    main()
