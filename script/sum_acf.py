#! /usr/bin/env python
from __future__ import print_function
import sys
import numpy

from cal_tc import WriterBase, ACFWriter, get_stringnames

import netCDF4 as netcdf

class TCSWriter(WriterBase):
    name = 'tcs'
    unit = '(kcal/mol)^2/fs'

def parse_options():
    """Parse and get the command line options."""

    # make argument parser
    import argparse
    parser = argparse.ArgumentParser(description=(
        'Average auto-correlation function over the given trajectories'))

    # add argument definitions
    parser.add_argument(dest='acf_fps', metavar='ACF_FILE',
            nargs='+',
            help=('The filepath of auto-correlation function data.'))

    parser.add_argument('-a', '--acf-file', metavar='FILE',
            dest='output_acf_fp', required=True,
            default='',
            help=('The filepath of acf data.'))

    parser.add_argument('-t', '--tcs-file', metavar='FILE',
            dest='tcs_fp', required=False,
            default='',
            help=('The filepath of tc time series data.'))

    parser.add_argument('-c', '--coefficient', metavar='COEFFICIENT',
            dest='coef', required=False,
            default=1.0, type=float,
            help='Multiply acf by given coefficient.')

    # make arguments
    return parser.parse_args()

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

def main():

    # get arguments
    opts = parse_options()

    import time
    t0 = time.time()

    # first step
    print('loading for {:>5}, '.format(1), end='')
    times, donors, acceptors, acf_sums = load_acf_first(opts.acf_fps[0])
    nacf = len(opts.acf_fps)
    nframe = len(times)
    dt = times[1] - times[0]

    t1 = time.time()
    print('time : {:>7.3f}(s)'.format(t1-t0))
    sys.stdout.flush()
    t0 = t1
    # average acf
    for iacf, acf_fp in enumerate(opts.acf_fps[1:],2):
        print('loading for {:>5}, '.format(iacf), end='')
        acfs = load_acf(acf_fp)
        acf_sums += acfs
        t1 = time.time()
        print('time : {:>7.3f}(s)'.format(t1-t0))
        sys.stdout.flush()
        t0 = t1

    acf_sums = acf_sums/float(nacf)

    # write averaged acf
    writer = ACFWriter(opts.output_acf_fp, nframe, dt)
    writer.write_all(donors, acceptors, acf_sums)

    # write time series transport coefficient
    if opts.tcs_fp:
        tcs_writer = TCSWriter(opts.tcs_fp, nframe, dt)
        tcs = cal_tcs(acf_sums, dt, opts.coef)
        tcs_writer.write_all(donors, acceptors, tcs)


if __name__ == '__main__':
    main()
