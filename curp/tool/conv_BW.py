#! /usr/bin/env python2
from __future__ import print_function
import sys, os

class BWDataParser:

    def __init__(self, filename):
        self.parse(filename)

    def parse(self, filename):
        file = open(filename, 'rb')
        self.rid_to_name_bw_pair = {}

        lines = ( line.strip() for line in file
                 if line.strip()
                 if not line.strip().startswith('%') 
                 if not line.strip().startswith('#') )


        for line in lines:

            rname, rid, bw_id = line.split()
            rid = int(rid)
            self.rid_to_name_bw_pair[rid] = (rname, bw_id)

        file.close()

if __name__ == '__main__':

    ec_file  = sys.stdin
    fn = sys.argv[1]

    parser = BWDataParser(fn)
    rid_to_name_bw_pair = parser.rid_to_name_bw_pair

    for line in ec_file:
        line = line.strip()
        if not line: continue

        if line.startswith('#'):
            continue

        cols = line.split()
        don_line, acc_line, rest = cols[0], cols[1], cols[2:]

        rid1 = don_line.split('_')[0]
        rid2 = acc_line.split('_')[0]

        if int(rid1) in rid_to_name_bw_pair:
            name1, bwid1 = rid_to_name_bw_pair.get(int(rid1))
            don_line = '{}_{}{} '.format(rid1, name1, bwid1)

        if int(rid2) in rid_to_name_bw_pair:
            name2, bwid2 = rid_to_name_bw_pair.get(int(rid2))
            acc_line = '{}_{}{} '.format(rid2, name2, bwid2)

        print(don_line, acc_line, '  '.join(rest))
