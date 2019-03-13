#! /usr/bin/env python2
from __future__ import print_function

import os, sys

def is_cap(rname):
    return rname.lower() in ['nme', 'ace']

if __name__ == '__main__':

    ec_file = sys.stdin

    lines = (line for line in ec_file
            if not line.startswith('#')
            if not line.isspace() )

    for line in lines:

        cols = line.split()
        don_line, acc_line, rest = cols[0], cols[1], cols[2:]

        rname1 = don_line.split('_')[1]
        rname2 = acc_line.split('_')[1]

        if is_cap(rname1): continue
        if is_cap(rname2): continue

        print(line.strip())

