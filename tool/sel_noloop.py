#! /usr/bin/env python2
from __future__ import print_function

import os, sys

def parse_ss(filename):
    """
    Parse the secondary structure file that is as follows:

    %TM1  1 25
    %TM2  40 90
    %SHEET1 102 110 
    .
    .
    .
    """

    file = open(filename, 'rb')
    lines = (line for line in file
            if line.startswith("%H")
            or line.startswith("%TM")
            or line.startswith("%S"))

            # if line.startswith("%TM")
            # if line.startswith("%H") )

    for line in lines:
        cols = line.split()
        name, ibeg, iend = cols[0].replace('%',''), cols[1], cols[2]
        yield int(ibeg), int(iend)

    file.close()

def is_loop(rid, rname, ibeg_iend_pairs, excluded_list):
    for ibeg, iend in ibeg_iend_pairs:
        if ibeg <= rid <= iend:
            return False
    else:
        return not include_resname(rname, excluded_list)

def include_resname(resname, rname_list):
    rname = resname.lower()
    for cs in rname_list:
        if rname.startswith(cs.lower()):
            return True
    else:
        return False
        
if __name__ == '__main__':

    ec_file = sys.stdin
    ss_fn = sys.argv[1]

    if len(sys.argv) >= 2:
        excluded_list = sys.argv[1:]
    else:
        excluded_list = ['WAT']

    ibeg_iend_pairs = list(parse_ss(ss_fn))

    lines = (line for line in ec_file
            if not line.startswith('#')
            if not line.isspace() )

    for line in lines:

        cols = line.split()
        don_line, acc_line, rest = cols[0], cols[1], cols[2:]

        rid1   = int(don_line.split('_')[0])
        rid2   = int(acc_line.split('_')[0])
        rname1 = don_line.split('_')[1]
        rname2 = acc_line.split('_')[1]

        if is_loop(rid1, rname1, ibeg_iend_pairs, excluded_list) \
                or is_loop(rid2, rname2, ibeg_iend_pairs, excluded_list):
            continue

        else:
            print(line.strip())
