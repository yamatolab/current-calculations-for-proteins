#! /usr/bin/env python2
import sys, os

if __name__ == '__main__':

    ec_file = sys.stdin

    if len(sys.argv) >= 2:
        excluded_list = sys.argv[1:]
    else:
        excluded_list = ['WAT']

    lines = (line for line in ec_file
            if not line.startswith("#")
            if not line.isspace() )

    for line in ec_file:

        cols = line.split()
        don_line, acc_line, rest = cols[0], cols[1], cols[2:]

        rid1 = int(don_line.split('_')[0])
        rid2 = int(acc_line.split('_')[0])
        rname1 = don_line.split('_')[1]
        rname2 = acc_line.split('_')[1]

        if rid1+1 != rid2: continue
        if rname1 in excluded_list or rname2 in excluded_list: continue

        print(line.strip())
