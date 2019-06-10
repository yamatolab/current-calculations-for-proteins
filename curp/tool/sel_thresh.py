#! /usr/bin/env python2
import sys, os

if __name__ == '__main__':

    # arguments
    if len(sys.argv) == 2:
        threshold   = float(sys.argv[1])
    else:
        print('{} < <ec_file> <threshold>'.format(sys.argv[0]))
        exit()

    ec_file = sys.stdin

    # remove the only space line and the comment line
    line_iter = ( line.strip() for line in ec_file if not line.isspace() )

    # store
    rkey_to_count = {}
    for line in line_iter:
        if line.startswith('#'):
            print(line)
            continue

        ec = float(line.split()[-1])
        if ec >= threshold: print(line)
