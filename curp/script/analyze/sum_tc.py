#! /usr/bin/env python
from __future__ import print_function
import sys
import numpy

class InvalidNumberPairs(Exception): pass

def write_tc(tcs, pairs, tcs_rms):
    """Write transport coefficient and dispersion to standard output."""
    for pair, tc, rms in zip(pairs, tcs, tcs_rms):
        don, acc = pair
        print('{:>12} {:>12}  {}  {}'.format(don, acc, tc, rms))

def load_tc(fp):
    fd = open(fp, 'rb')
    lines = (line for line in fd
            if not line.startswith('#')
            if not line.isspace())

    pairs = []
    ipair_to_tc = []
    for line in lines:
        cols = line.split()
        don, acc = cols[0], cols[1]
        tc = float(cols[2])

        pairs += [(don, acc)]
        ipair_to_tc += [tc]

    return pairs, ipair_to_tc

def summarize_tc(tc_fps, **kwds):
    """Summarize transport coefficient"""

    ipair_to_tc_ary = []
    npair = None
    for fp in tc_fps:
        pairs, ipair_to_tc = load_tc(fp)
        ipair_to_tc_ary += [ipair_to_tc]

        if npair is None: npair = len(pairs)

        if npair != len(pairs):
            msg = "Number of pairs must be {}, but got {} in {}"
            raise InvalidNumberPairs(msg
                    .format(npair, len(pairs), fp))

    tc_ary = numpy.array(ipair_to_tc_ary)

    tc_avg = numpy.average(tc_ary, 0)
    tc_rms = numpy.std(tc_ary, 0)

    write_tc(tc_avg, pairs, tc_rms)

if __name__ == '__main__':
    from curp.console import arg_sum_tc, exec_command

    parser = arg_sum_tc()
    exec_command(parser)
