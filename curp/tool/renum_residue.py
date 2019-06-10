#! /usr/bin/env python2
from __future__ import print_function

import sys, os


def parse_line(line):
    cols = line.split()
    don_line, acc_line, rest = cols[0], cols[1], cols[2:]

    rid_don, rname_don = don_line.split('_')
    rid_don = int(rid_don)
    rid_acc, rname_acc = acc_line.split('_')
    rid_acc = int(rid_acc)

    return rid_don, rid_acc, rname_don, rname_acc, rest

if __name__ == '__main__':

    # arguments
    if len(sys.argv) >= 2:
        renums = sys.argv[1:]
    else:
        print('{} < <ec_file> <1:15> <250:300> ...'.format(sys.argv[0]))
        exit()


    # remove the only space line and the comment line
    ec_file = sys.stdin
    lines = [ line.strip() for line in ec_file
            if not line.isspace() if not line.startswith('#') ]

    # obtain maximum residue number
    rid_max = 0
    for line in lines:
        rid1, rid2, rname1, rname2, rest = parse_line(line)
        rid_max = max(rid_max, rid1, rid2)

    new_rids = range(rid_max+1)
    old_rids = range(rid_max+1)

    # make dictionary to alter the given residue identifiers
    for renum in renums:
        rid_old, rid_new  = renum.split(':')
        rid_old, rid_new = int(rid_old), int(rid_new)

        shift = rid_new - rid_old
        new_rids[rid_old:] = [ rid+shift for rid in old_rids[rid_old:] ]

    # convert residue identifiers into new ones
    for line in lines:
        rid1, rid2, rname1, rname2, rest = parse_line(line)
        ec = rest[0]

        # print
        fmt = '{:05}_{:<7}  '
        print(fmt.format(new_rids[rid1], rname1),
                fmt.format(new_rids[rid2], rname2), ec)
