#! /usr/bin/env python2

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
    if len(sys.argv) == 2:
        sort_method = sys.argv[1]
    else:
        print('{} < <ec_file> <sort method>'.format(sys.argv[0]))
        exit()

    ec_file = sys.stdin

    # remove the only space line and the comment line
    line_iter = ( line.strip() for line in ec_file
            if not line.isspace() if not line.startswith('#') )

    # store
    rkey_to_count = {}
    for line in line_iter:
        rid1, rid2, rname1, rname2, rest = parse_line(line)

        ec = float(rest[0])

        count1 = rkey_to_count.get((rid1, rname1))
        rkey_to_count[rid1, rname1] = count1+1 if count1 else 1

        count2 = rkey_to_count.get((rid2, rname2))
        rkey_to_count[rid2, rname2] = count2+1 if count2 else 1

    # sort
    # seq, count, rev_count
    if sort_method == 'seq':
        sorted_items = sorted( rkey_to_count.items(),
            key=lambda rkey: rkey[0][0], reverse=False)
    elif sort_method == 'count':
        sorted_items = sorted( rkey_to_count.items(),
            key=lambda rkey: rkey[1], reverse=True)
    elif sort_method == 'rev_count':
        sorted_items = sorted( rkey_to_count.items(),
            key=lambda rkey: rkey[1], reverse=False)
    else:
        sorted_items = rkey_to_count.items()

    # print
    print('#{:>8}{:>5}   {:>2}'.format('resname', 'rid', 'count'))
    for (rid, rname), count in sorted_items:
        print(' {:>8}{:>5}   {:>2}'.format(rname, rid, count))

