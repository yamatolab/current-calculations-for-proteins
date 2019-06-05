#! /usr/bin/env python
from __future__ import print_function
import sys
import numpy

from cal_tc import WriterBase, ACFWriter, get_stringnames

import netCDF4 as netcdf

class TCSWriter(WriterBase):
    name = 'tcs'
    unit = '(kcal/mol)^2/fs'

def load_acf(acf_fp, dataname='acf'):

    ncfile  = netcdf.Dataset(acf_fp, mode='r')
    acfs  = ncfile.variables[dataname][:]
    ncfile.close()

    return acfs

def load_acf_first(acf_fp, dataname='acf'):

    ncfile  = netcdf.Dataset(acf_fp, mode='r')

    donors    = get_stringnames( ncfile.variables['donors'][:] )
    acceptors = get_stringnames( ncfile.variables['acceptors'][:] )

    times     = ncfile.variables['time'][:]
    acfs       = ncfile.variables[dataname][:]

    ncfile.close()

    return times, donors, acceptors, acfs

def cal_tcs(acf, dt, coef=1.0):
    """Calculate time series of transport coefficient."""
    return coef * dt*1000.0 * numpy.add.accumulate(acf, 1)

def sum_acf(acf_fps, output_acf_fp, tcs_fp='', coef=1., **kwds):
    import time
    t0 = time.time()

    # first step
    print('loading for {:>5}, '.format(1), end='')
    times, donors, acceptors, acf_sums = load_acf_first(acf_fps[0])
    nacf = len(acf_fps)
    nframe = len(times)
    dt = times[1] - times[0]

    t1 = time.time()
    print('time : {:>7.3f}(s)'.format(t1-t0))
    sys.stdout.flush()
    t0 = t1
    # average acf
    for iacf, acf_fp in enumerate(acf_fps[1:],2):
        print('loading for {:>5}, '.format(iacf), end='')
        acfs = load_acf(acf_fp)
        acf_sums += acfs
        t1 = time.time()
        print('time : {:>7.3f}(s)'.format(t1-t0))
        sys.stdout.flush()
        t0 = t1

    acf_sums = acf_sums/float(nacf)

    # write averaged acf
    writer = ACFWriter(output_acf_fp, nframe, dt)
    writer.write_all(donors, acceptors, acf_sums)

    # write time series transport coefficient
    if tcs_fp:
        tcs_writer = TCSWriter(tcs_fp, nframe, dt)
        tcs = cal_tcs(acf_sums, dt, coef)
        tcs_writer.write_all(donors, acceptors, tcs)


if __name__ == '__main__':
    main()
