#! /usr/bin/env python2
from __future__ import print_function

import os, sys

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

        
if __name__ == '__main__':

    ec_file = sys.stdin
    helix_fn = sys.argv[1]
    is_lt_4 = len(sys.argv)==3 and sys.argv[2]=='1'

    ibeg_iend_pairs = list(parse_helix(helix_fn))

    for line in ec_file:
        if line.startswith('#'):
            print(line.strip())
            continue

        cols = line.split()
        don_line, acc_line, rest = cols[0], cols[1], cols[2:]

        rid1 = int(don_line.split('_')[0])
        rid2 = int(acc_line.split('_')[0])

        # if rid1+4 == rid2:
        if is_lt_4:
            if rid1+4 >= rid2:
                for ibeg, iend in ibeg_iend_pairs:

                    if ibeg<=rid1<=iend and ibeg<=rid2<=iend:
                        # print(rid1, rid2, ibeg, iend)
                        break

                else:
                    print(line.strip())

            else:
                print(line.strip())

        else:
            if rid1+4 == rid2:
                for ibeg, iend in ibeg_iend_pairs:

                    if ibeg<=rid1<=iend and ibeg<=rid2<=iend:
                        # print('this', rid1, rid2, ibeg, iend)
                        break

                else:
                    print(line.strip())

            else:
                print(line.strip())


