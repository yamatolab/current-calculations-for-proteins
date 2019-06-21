#! /usr/bin/env python2
from __future__ import print_function

import os, sys

def parse_asa(filename):
    """
    1
    2
    3
    5
    9
    .
    .
    341
    342
    """

    file = open(filename, 'rb')

    for line in file:
        if line.startswith('#'): continue
        yield int(line.strip())

    file.close()


def main(ec_fn, asa_fn, **kwds):
    if ec_fn is None:
        ec_file = sys.stdin
    else:
        ec_file = open(ec_fn, 'r')
        lines = (line for line in ec_file
                if not line.startswith('#')
                if not line.isspace() )
    ec_file.close()

    for line in lines:
        # if line.startswith('#'):
            # print(line.strip())
            # continue
        # print(line)

        cols = line.split()
        don_line, acc_line, rest = cols[0], cols[1], cols[2:]

        rid1 = int(don_line.split('_')[0])
        rid2 = int(acc_line.split('_')[0])

        if rid1 in asa_resids or rid2 in asa_resids:
        # if rid1 in asa_resids and rid2 in asa_resids:
            # print(rid1, rid2)
            continue

        else:
            print(line.strip())

if __name__ == '__main__':
    from curp.tool.sele.console import exec_command, arg_sel_noasa
    parser = arg_sel_noasa()
    exec_command(parser)
