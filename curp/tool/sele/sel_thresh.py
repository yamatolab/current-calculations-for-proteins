#! /usr/bin/env python2
import sys
import os


def main(ec_fn, threshold, **kwds):
    # remove the only space line and the comment line
    if ec_fn is None:
        ec_file = sys.stdin
    else:
        ec_file = open(ec_fn, 'r')
    line_iter = (line.strip() for line in ec_file
                 if not line.isspace())
    ec_file.close()

    # store
    rkey_to_count = {}
    for line in line_iter:
        if line.startswith('#'):
            print(line)
            continue

        ec = float(line.split()[-1])
        if ec >= threshold: print(line)


if __name__ == '__main__':
    from curp.tool.sele.console import exec_command, arg_sel_thresh
    parser = arg_sel_thresh()
    exec_command(parser)
