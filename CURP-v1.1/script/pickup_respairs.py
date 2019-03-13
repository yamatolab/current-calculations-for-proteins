#! /usr/bin/env python
from __future__ import print_function

import os, sys
import numpy
import argparse

import curp_module

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
    tpl = curp_module.get_topology(args.prmtop_fmt, args.prmtop_fn)
    natom = tpl.get_natom()

    # generate residue group
    res_info  = tpl.get_residue_info()
    atom_info = tpl.get_atom_info()
    rname_iatoms_pairs = list( (rname, list(numpy.array(iatoms)))
            for rname, iatoms in gen_residue_group(
                res_info, atom_info, name_fmt=args.format))

    resnames = [ resname for resname, iatoms in rname_iatoms_pairs ]
    nres = len(resnames)

    # get trajectory
    crd_parser = curp_module.gen_crds(args.trj_fns, args.input_trj_fmt,
            natom=natom, fst_lst_int=(1,-1,args.interval), use_pbc=False)
    crds = ( crd for itraj, crd, box in crd_parser )

    # generate residue pair table
    import lib_pickup
    # cutoff_method = dict(
            # com      = lib_pickup.is_com,
            # nearest  = lib_pickup.is_nearest,
            # farthest = lib_pickup.is_farthest)[args.cutoff_method]


    res_beg_end_pairs = [ (min(iatoms), max(iatoms)) 
            for rname, iatoms in rname_iatoms_pairs ]

    lib_pickup.pickup.res_beg_end_pairs = numpy.array(res_beg_end_pairs)
    lib_pickup.pickup.cutoff2 = cutoff2
    if args.cutoff_method == 'com':
        lib_pickup.pickup.cutoff_method = 1
    elif args.cutoff_method == 'nearest':
        lib_pickup.pickup.cutoff_method = 2
    elif args.cutoff_method == 'farthest':
        lib_pickup.pickup.cutoff_method = 3
    else:
        pass

    respair_table_traj = (gen_respair_table_fortran(crd, nres, resnames)
            for crd in crds )

    respair_iter = gen_respair_table_over_traj(
            respair_table_traj, args.is_union)

    # setup excluded resids
    def gen_excluded_respair(respairs, excluded_resids):
        ext_resids = []
        for ext_resid in excluded_resids:
            rid_beg, rid_end = ext_resid.split(':')
            rid_beg, rid_end = int(rid_beg), int(rid_end)
            ext_resids.extend( range(rid_beg, rid_end+1) )

        # print(ext_resids)

        flag_i = False
        for rname_i, rnames_j in respairs:
            rid_i = int(rname_i.split('_')[0])
            if rid_i in ext_resids: flag_i = True

            new_rnames_j = []

            for rname_j in rnames_j:
                rid_j = int(rname_j.split('_')[0])
                if flag_i and rid_j in ext_resids: continue

                new_rnames_j.append(rname_j)

            else:
                flag_i = False

            if new_rnames_j:
                yield rname_i, new_rnames_j

    def gen_excluded_respair_both(respairs, excluded_resids):
        ext_resids = []
        for ext_resid in excluded_resids:
            rid_beg, rid_end = ext_resid.split(':')
            rid_beg, rid_end = int(rid_beg), int(rid_end)
            ext_resids.extend( range(rid_beg, rid_end+1) )

        # print(ext_resids)

        for rname_i, rnames_j in respairs:
            rid_i = int(rname_i.split('_')[0])
            if rid_i in ext_resids: continue

            new_rnames_j = []

            for rname_j in rnames_j:
                rid_j = int(rname_j.split('_')[0])
                if rid_j in ext_resids: continue

                new_rnames_j.append(rname_j)

            if new_rnames_j:
                yield rname_i, new_rnames_j

    # write residue pair data
    print('# residue pair within {:4.2f}'.format(cutoff))
    npair = 0
    if args.ext_resids and not args.is_both:

        for rname_i, rnames_j in gen_excluded_respair(
                respair_iter, args.ext_resids):

            if rnames_j:
                npair += len(rnames_j)
                print('[{}]'.format(rname_i))
                print('  '.join(rnames_j))

    elif args.ext_resids and args.is_both:

        for rname_i, rnames_j in gen_excluded_respair_both(
                respair_iter, args.ext_resids):

            if rnames_j:
                npair += len(rnames_j)
                print('[{}]'.format(rname_i))
                print('  '.join(rnames_j))

    else:
        for rname_i, rnames_j in respair_iter:

            if rnames_j:
                npair += len(rnames_j)
                print('[{}]'.format(rname_i))
                print('  '.join(rnames_j))


    print('# number of pairs = ', npair)

def get_arguments():

    parser = argparse.ArgumentParser(description=
            'Pick up the residue pair table to given file as a ndx format.')

    # add argument definitions
    parser.add_argument(
            dest='trj_fns', nargs='+',
            help='specify the trajectory files.')

    parser.add_argument('-if', '--input-format', metavar='FORMAT',
            dest='input_trj_fmt', required=True,
            help='specify the format of the trajectory.')

    parser.add_argument(
            '-p', '--input-prmtop-file', dest='prmtop_fn', required=True,
            help='specify the prmtop file to be amber format.')

    parser.add_argument(
            '-pf', '--input-prmtop-format', dest='prmtop_fmt', required=True,
            help='specify the format of the prmtop file.')
    
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

    parser.add_argument('-e', '--exclude-resids', metavar='FIRST:LAST',
            dest='ext_resids', required=False,
            default='', nargs='*',
            help=('residues that you want to exclude.'+
                'ex) 5:10 80:150' ))

    parser.add_argument(
            '-b', '--enable-residues-both-cxcluded', dest='is_both',
            required=False, action='store_true',
            help=('whether include residue pairs that exclude residues'))


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

def gen_respair_table_fortran(crd, nres, resnames):

    import lib_pickup
    cands = lib_pickup.pickup.get_candidate_pairs(crd, nres)

    for ires_1 in range(nres-1):
        jreslist = []
        for jres_1 in range(ires_1+1, nres):
            if cands[ires_1, jres_1]: jreslist.append( resnames[jres_1] )

        yield resnames[ires_1], jreslist

def gen_respair_table(crd, resnames, rname_iatoms_pairs, cutoff_method,cutoff2):

    res_crds = [ crd[iatoms] for rname, iatoms in rname_iatoms_pairs ]

    for i, (rname_i, rcrd_i) in enumerate(zip(resnames, res_crds)):
        resj_cands = []
        for j, (rname_j, rcrd_j) in enumerate(zip(resnames, res_crds)):
            if j <= i: continue

            if cutoff_method(rcrd_i, rcrd_j, cutoff2):
                resj_cands.append(rname_j)

        yield rname_i, resj_cands


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

    for itrj, respair_table in enumerate(respair_table_traj,1):
        print('# itraj = ', itrj)

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
