#! /usr/bin/env python
from __future__ import print_function

import numpy

from curp import get_tpl
from curp.script.conv_trj import gen_trj
from curp.script.analyze.lib_pickup import pickup
from curp.table.group import gen_residue_group

def pickup_respairs(trj_fns, input_trj_fmt, prmtop_fn, prmtop_fmt, interval=1,
         cutoff_method='nearest', cutoff=5., trim_resnames=[],
         format='{rid:05}_{rname}', is_union=True, ext_resids='',
         is_both=False, **kwds):
    """Pick up the residue pair table to given file as a ndx format."""

    # Get parameters
    cutoff = float(cutoff)
    cutoff2 = cutoff*cutoff

    # Get topology
    tpl = get_tpl(prmtop_fmt, prmtop_fn)
    natom = tpl.get_natom()

    # Generate residue group
    res_info  = tpl.get_residue_info()
    atom_info = tpl.get_atom_info()
    rname_iatoms_pairs = list( (rname, list(numpy.array(iatoms)))
            for rname, iatoms in gen_residue_group(
                res_info, atom_info, name_fmt=format))

    resnames = [ resname for resname, iatoms in rname_iatoms_pairs ]
    nres = len(resnames)

    # Get trajectory
    crd_parser = gen_trj(tpl, trj_fns, input_trj_fmt,
            trj_type='crd', input_fst_lst_int=[(1, -1, interval)],
            use_pbc=False)
    crds = ( crd for itraj, crd, box in crd_parser )

    res_beg_end_pairs = [ (min(iatoms), max(iatoms))
            for rname, iatoms in rname_iatoms_pairs ]

    pickup.res_beg_end_pairs = numpy.array(res_beg_end_pairs)
    pickup.cutoff2 = cutoff2
    if cutoff_method == 'com':
        pickup.cutoff_method = 1
    elif cutoff_method == 'nearest':
        pickup.cutoff_method = 2
    elif cutoff_method == 'farthest':
        pickup.cutoff_method = 3
    else:
        pass

    respair_table_traj = (gen_respair_table_fortran(crd, nres, resnames)
            for crd in crds )

    respair_iter = gen_respair_table_over_traj(
            respair_table_traj, is_union)

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
    if ext_resids and not is_both:

        for rname_i, rnames_j in gen_excluded_respair(
                respair_iter, ext_resids):

            if rnames_j:
                npair += len(rnames_j)
                print('[{}]'.format(rname_i))
                print('  '.join(rnames_j))

    elif ext_resids and is_both:

        for rname_i, rnames_j in gen_excluded_respair_both(
                respair_iter, ext_resids):

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


def cal_dist2(p1, p2):
    # return sum( (p1 - p2)**2 )
    return (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2

def zip_double(data_iter):

    for i, e1 in enumerate(data_iter):
        for j, e2 in enumerate(data_iter):
            if j <= i: continue
            yield e1, e2

def gen_respair_table_fortran(crd, nres, resnames):

    cands = pickup.get_candidate_pairs(crd, nres)

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
    from curp.console import arg_respairs, exec_command

    parser = arg_respairs()
    exec_command(parser)
