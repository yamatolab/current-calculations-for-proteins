from __future__ import print_function

import os, sys

import numpy
topdir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)
from exception import CurpException

class MaximumDepthError(CurpException): pass

def gen_searched_wells(rdfs, dr, increment=5, eps=0.0000001):
    for iatm_1, rdf in enumerate(rdfs):
        try:
            r_max,rdf_max,r_min,rdf_min = search_well(rdf, dr, increment, eps)
        except MaximumDepthError as e:
            msg = 'Reached maximum depth: {} at the atom index {}.'
            raise MaximumDepthError(msg.format(e.num_rdfs, iatm_1+1))

        yield r_max, rdf_max, r_min, rdf_min

def search_well(rdf, dr, increment=5, eps=0.0000001):

    # get a index that rdf is maximum.
    num_rdfs = len(rdf)
    rdf_max  = max(rdf)
    i_max = index_float_array(rdf, rdf_max)

    for i in range(num_rdfs):
        i_next = i_max+1 + (i+1)*increment
        rdf_next = rdf[i_max+1:i_next]

        # get minimum and its index
        rdf_min = min(rdf_next)
        i_min = index_float_array(rdf_next, rdf_min)
        i_min += i_max+1

        # judgement
        if num_rdfs == i_min+1:
            break
        else:
            if rdf_min <= rdf[i_min+1]:
            # if rdf_min == 0.0 or rdf_min <= rdf[i_min+1]:
                break

    else:
        error_obj = MaximumDepthError()
        error_obj.num_rdfs = num_rdfs
        raise error_obj

    return dr*(i_max+1), rdf_max, dr*(i_min+1), rdf_min

def index_float_array(array, value, eps=0.0000001):
    criterion = (value-eps <= array) & (array <= value+eps)
    index_item = numpy.where(criterion)
    if len(index_item[0]) == 0:
        return None
    else:
        return index_item[0][0]


if __name__ == '__main__':

    import os, sys
    topdir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '..', 'parser', 'amber'))
    if topdir not in sys.path:
        sys.path.insert(0, topdir)
    import trajectory
    from calrdf import average_rdf
    
    def print_wells(info):
        rmax, rdf_max, rmin, rdf_min = info
        msg = '{rmax:4.2f} {rdf_max:10.5f}   {rmin:4.2f} {rdf_min:10.5f}'
        print(msg.format(rmax=rmax, rdf_max=rdf_max, rmin=rmin,rdf_min=rdf_min))

    crd_fn = '../parser/amber/test/sam-nwat.mdcrd'
    parser = trajectory.CoordinateParser(crd_fn, 1799)

    rmax = 3.0
    dr   = 0.01

    from benchmarker import Benchmarker
    with Benchmarker(width=30) as bm:

        with bm('calculate average rdf'):
            rdfs = average_rdf(parser, rmax=rmax, dr=dr, interval=1,
                    average=True, per_area=True) 
        # with bm('printing average rdf'):
        #     for iatm_1, rdf in enumerate(rdfs):
        #         print()
        #         print(iatm_1+1)
        #         for rd in rdf:
        #             print(rd)

        with bm('search well of rdf'):
            for iatm_1, tup in enumerate(
                    gen_searched_wells(rdfs, dr, increment=5)):
                print(iatm_1+1)
                print_wells(tup)


