from __future__ import print_function

import math
import numpy

def cal_rdf(crd, num_rs, rmax=5.0, dr=0.1, per_area=True):
    """Calclate a radial distribution function on each atom."""
    natom = len(crd)
    rdfs = numpy.zeros((natom, num_rs)) # rdf on each atom.
    
    for rdf_i, r_i in zip(rdfs, crd):
        for r_j in crd:
            r_ij = r_i - r_j
            l_ij = math.sqrt(numpy.dot(r_ij, r_ij))
            if l_ij <= 0.1: continue
            if l_ij >= rmax: continue
            rindex = int(round(l_ij/dr)) - 1
            rdf_i[rindex] += 1

    if per_area:
        for rindex in range(num_rs):
            r = dr * (rindex + 1)
            rdfs[:, rindex] = rdfs[:, rindex] / (r*r)
                
    return rdfs

# replace the calrdf routine to fortran one
import lib_calrdf
cal_rdf = lib_calrdf.calrdf

def average_rdf(parser, rmax=5.0, dr=0.1
        , interval=1, average=True, per_area=True):
    """Calculate a average rdf on each atom."""
    num_rs = int(math.ceil(rmax/dr))
    ntraj = 0

    # calculate first rdfs
    i = 0
    for istep, (crd, box) in parser:
        i += 1
        istep, (crd, box) = parser.next()

        if i == interval:
            ntraj += 1
            break

    natom = len(crd)
    rdfs_total = cal_rdf(crd, num_rs, rmax=rmax, dr=dr, per_area=False)

    # calclate rest rdfs
    i = 0
    for istep, (crd, box) in parser:
        i += 1

        if i == interval:
            ntraj += 1
            rdfs = cal_rdf(crd, num_rs, rmax=rmax, dr=dr, per_area=False)
            rdfs_total += rdfs
            i = 0

    if per_area:
        for rindex in range(num_rs):
            r = dr * (rindex + 1)
            rdfs_total[:, rindex] = rdfs_total[:, rindex] / (r*r)

    if average:
        return rdfs_total/ntraj
    else:
        return rdfs_total

if __name__ == '__main__':
    import os, sys
    topdir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '..', 'parser', 'amber'))
    if topdir not in sys.path:
        sys.path.insert(0, topdir)
    import trajectory

    crd_fn = '../parser/amber/test/sam-nwat.mdcrd'
    parser = trajectory.CoordinateParser(crd_fn, 1799)

    rmax = 2.5
    dr   = 0.01

    from benchmarker import Benchmarker
    with Benchmarker(width=30) as bm:

        # with bm('python code: calculation'):

        #     rdfs_py = cal_rdf(crd, num_rs, rmax=rmax, dr=dr, per_area=True)

        # with bm('fortran code: calculation'):

        #     rdfs_for = lib_calrdf.calrdf(
        #             crd, num_rs, rmax=rmax, dr=dr, per_area=True)

        # with bm('printing'):
        #     for iatm_1, (rdf_py, rdf_for) in enumerate(zip(rdfs_py, rdfs_for)):
        #         print()
        #         print(iatm_1+1, rdf_py-rdf_for)

        with bm('calculate average rdf'):
            rdfs = average_rdf(parser, rmax=rmax, dr=dr, interval=1,
                    average=False, per_area=False) 
        with bm('printing average rdf'):
            for iatm_1, rdf in enumerate(rdfs):
                print()
                print(iatm_1+1)
                for rd in rdf:
                    print(rd)

