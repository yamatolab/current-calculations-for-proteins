#! /usr/bin/env python2
from __future__ import print_function
import os
import sys

def parse_helix(filename):
    """
    %TM1  1 25
    %TM2  40 90
    .
    .
    .
    """

    file = open(filename, 'rb')
    lines = (line for line in file
            if line.startswith("%H")
            or line.startswith("%TM") )

            # if line.startswith("%TM")
            # if line.startswith("%H") )

    for line in lines:
        cols = line.split()
        name, ibeg, iend = cols[0].replace('%',''), cols[1], cols[2]
        yield int(ibeg), int(iend)

    file.close()

def main(ec_fn, helix_fn, arg_3=None, **kwds):

    ibeg_iend_pairs = list(parse_helix(helix_fn))

    if ec_fn is None:
        ec_file = sys.stdin
    else:
        ec_file = open(ec_fn, 'r')

    for line in ec_file:
        if line.startswith('#'):
            print(line.strip())
            continue

        cols = line.split()
        don_line, acc_line, rest = cols[0], cols[1], cols[2:]
        rid1 = int(don_line.split('_')[0])
        rid2 = int(acc_line.split('_')[0])

        if arg_3 == '1':
            if rid1+4 >= rid2:
                for ibeg, iend in ibeg_iend_pairs:
                    if ibeg<=rid1<=iend and ibeg<=rid2<=iend:
                        print(line.strip())

        else:
            if rid1+4 == rid2:
                for ibeg, iend in ibeg_iend_pairs:
                    if ibeg<=rid1<=iend and ibeg<=rid2<=iend:
                        print(line.strip())

    ec_file.close()


if __name__ == '__main__':
    from curp.tool.sele.console import exec_command, arg_sel_helix
    parser = arg_sel_helix()
    exec_command(parser)
