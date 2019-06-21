from __future__ import print_function
import sys
import numpy

from curp.script.cal_tc import get_stringnames
from curp.script.sum_acf import load_acf_first

import netCDF4 as netcdf



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


def get_ncdata(acf_fp, group_ranges, dataname='acf', prefix='acf', **kwds):
    """Get simple text data from file in netcdf format by given arguments."""

    index_iter = gen_indices(group_ranges)

    # first step
    times, donors, acceptors, acfs = load_acf_first(acf_fp, dataname)

    for i in index_iter:
        i_1 = i-1
        donor, acceptor, acf = donors[i_1], acceptors[i_1], acfs[i_1]
        write_data(prefix, donor, acceptor, acf, times)


if __name__ == '__main__':
    from curp.console import arg_ncdata, exec_command

    parser = arg_ncdata()
    exec_command(parser)
