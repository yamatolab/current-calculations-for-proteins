#! /usr/bin/env python2
import sys
import os


def main(ec_fn=None, excluded_list=['WAT'], **kwds):

    if ec_fn is None:
        ec_file = sys.stdin
    else:
        ec_file = open(ec_fn, 'r')
    lines = (line for line in ec_file
             if not line.startswith("#")
             if not line.isspace())

    for line in lines:
        cols = line.split()
        don_line, acc_line, rest = cols[0], cols[1], cols[2:]

        rid1 = int(don_line.split('_')[0])
        rid2 = int(acc_line.split('_')[0])
        rname1 = don_line.split('_')[1]
        rname2 = acc_line.split('_')[1]

        if rid1+1 != rid2: continue
        if rname1 in excluded_list or rname2 in excluded_list: continue

        print(line.strip())


if __name__ == '__main__':
    from curp.tool.sele.console import exec_command, arg_sel_neighbor
    parser = arg_sel_neighbor()
    exec_command(parser)
