#! /usr/bin/env python
from __future__ import print_function
import sys
import numpy

from cal_tc import get_stringnames
from sum_acf import load_acf_first

import netCDF4 as netcdf

def parse_options():
    """Parse and get the command line options."""

    # make argument parser
    import argparse
    parser = argparse.ArgumentParser(description=(
        'Get simple text data from file in netcdf format by given arguments.'))

    # add argument definitions
    parser.add_argument(dest='acf_fp', metavar='ACF_FILE',
            help=('The filepath of auto-correlation function data.'))

    parser.add_argument('-r', '--group-ranges', metavar='FRIST:LAST',
            nargs='*', 
            dest='group_ranges', required=True,
            default='',
            help='The pair range list to want to gain.')

    parser.add_argument('-n', '--dataname',
            dest='dataname', required=False,
            default='acf',
            help='The name of the netcdf data you want to gain.')

    parser.add_argument('-o', '--output-prefix',
            dest='prefix', required=False,
            default='acf',
            help=('The prefix of the files to want to write, '
                +'that includes directory path.'))

    # make arguments
    return parser.parse_args()

def gen_indices(grange_lines):
    for grange_line in grange_lines:
        fst, lst = grange_line.split(':')
        if float(int(fst))<=0:
            msg = 'The start index is invalid: {}.'.format(fst)
            raise Exception(msg)

        for index in range(int(fst), int(lst)+1):
            yield index

def write_data(prefix, donor, acceptor, acf, times):
    fn = '{:}-{:}-{:}.dat'.format(prefix, donor, acceptor)
    numpy.savetxt(fn, numpy.array([times, acf]).T, fmt='%12.4f   % 12.7e')

def main():

    # get arguments
    opts = parse_options()

    index_iter = gen_indices(opts.group_ranges)

    # first step
    times, donors, acceptors, acfs = load_acf_first(opts.acf_fp, opts.dataname)

    for i in index_iter:
        i_1 = i-1
        donor, acceptor, acf = donors[i_1], acceptors[i_1], acfs[i_1]
        write_data(opts.prefix, donor, acceptor, acf, times)

if __name__ == '__main__':
    main()
