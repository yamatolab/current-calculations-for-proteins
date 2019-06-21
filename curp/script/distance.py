#! /usr/bin/env python
from __future__ import print_function

import numpy

from curp.script.lib_distance import distance

def do_distance(tpl, trj, trj_type,
                dist_format=None, dist_cutoff=5.0, dist_method='com', **kwds
                ):

    if not trj_type == 'crd':
        msg = 'Distance method is valid only for coordinate trajectory.'
        raise Exception(msg)

    from curp.table.group import gen_residue_group

    # Generate residue group
    res_info  = tpl.get_residue_info()
    atom_info = tpl.get_atom_info()
    rname_iatoms_pairs = list( (rname, list(numpy.array(iatoms)))
            for rname, iatoms in gen_residue_group(
                res_info, atom_info, name_fmt=dist_format))

    resnames = [ resname for resname, iatoms in rname_iatoms_pairs ]
    nres = len(resnames)

    res_beg_end_pairs = [ (min(iatoms), max(iatoms))
            for rname, iatoms in rname_iatoms_pairs ]

    mod = distance
    mod.res_beg_end_pairs = numpy.array(res_beg_end_pairs)
    mod.cutoff2 = dist_cutoff ** 2
    if dist_method == 'com':        mod.method = 1
    elif dist_method == 'nearest':  mod.method = 2
    elif dist_method == 'farthest': mod.method = 3
    else: pass

    dist2s_trj = (get_dist2s(crd, nres) for ifrm, crd, box in trj )

    # Calculate average and rms of the distances
    dists_min, dists_avg, dists_max, dists_rms = get_min_avg_max_rms(dist2s_trj)

    # Write data
    write_distance(dists_min, dists_avg, dists_max, dists_rms, resnames)

def get_dist2s(crd, nres):
    dist2s = distance.cal_dist2s(crd, nres)
    return dist2s

def get_min_avg_max_rms(dist2s_trj):

    import time
    t0 = time.time()
    print('# itraj = {:>8}'.format(1), end='')
    dist2s_fst = dist2s_trj.next()
    dist2s_sum = dist2s_fst
    dists_sum  = numpy.sqrt(dist2s_fst)

    dist2s_min = dist2s_fst.copy()
    dist2s_max = dist2s_fst.copy()

    t1 = time.time()
    print('{:12.4f}'.format(t1-t0))
    t0 = t1
    for itrj, dist2s in enumerate(dist2s_trj, 2):
        print('# itraj = {:>8}'.format(itrj), end='')

        dist2s_sum += dist2s
        dists_sum  += numpy.sqrt(dist2s)
        dist2s_min = numpy.minimum(dist2s_min, dist2s)
        dist2s_max = numpy.maximum(dist2s_max, dist2s)

        t1 = time.time()
        print('{:12.4f}'.format(t1-t0))
        t0 = t1

    nframe = itrj

    dists_avg = dists_sum/float(nframe)
    dists_rms = numpy.sqrt(dist2s_sum/float(nframe) - dists_avg**2)
    return numpy.sqrt(dist2s_min), dists_avg, numpy.sqrt(dist2s_max), dists_rms

def write_distance(dmin, davg, dmax, drms, resnames):

    nres = len(resnames)
    '# write distance data'
    header_fmt = '# {:>12} {:>12} {:>12s} {:>12s} {:>12s} {:>12s}'
    data_fmt   = '  {:>12} {:>12} {:>12.3f} {:>12.3f} {:>12.3f} {:>12.3f}'
    print(header_fmt.format(
        'residue 1', 'residue 2', 'min_dist', 'avg_dist', 'max_dist', 'rms'))
    for ires_1 in range(nres-1):
        res_i = resnames[ires_1]
        for jres_1 in range(ires_1+1, nres):
            res_j = resnames[jres_1]
            print(data_fmt.format(
                res_i, res_j, dmin[ires_1,jres_1], davg[ires_1,jres_1],
                dmax[ires_1,jres_1], drms[ires_1,jres_1]))

    print('# Finished completely.')

def load_mask(fn):

    with open(fn, 'rb') as file:
        ids_str = (line.split()[1] for line in file
                if line.startswith('ATOM'))

        ids = [ int(id_str)-1 for id_str in ids_str ]

    return ids

if __name__ == '__main__':
    pass
