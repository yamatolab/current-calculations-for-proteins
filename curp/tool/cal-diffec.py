#! /usr/bin/env python2
from __future__ import print_function

import sys, os

def main():

    opts = parse_options()

    tc_fp1   = opts.tc_fps[0]
    tc_fp2   = opts.tc_fps[1]
    tc_outfp = opts.output_fp

    pair_to_tc1 = { (don,acc):tc for don,acc,tc,rest in load_tcdata(tc_fp1) }
    pair_to_tc2 = { (don,acc):tc for don,acc,tc,rest in load_tcdata(tc_fp2) }

    ipair_to_pair1 = { (get_id(don),get_id(acc)):(don, acc)
            for (don,acc) in pair_to_tc1.keys() }

    ipair_to_pair2 = { (get_id(don),get_id(acc)):(don, acc)
            for (don,acc) in pair_to_tc2.keys() }

    ipair_to_tc1 = { (get_id(don),get_id(acc)):tc
            for (don,acc), tc in pair_to_tc1.items() }
    ipair_to_tc2 = { (get_id(don),get_id(acc)):tc
            for (don,acc), tc in pair_to_tc2.items() }

    ipair1 = set(ipair_to_pair1)
    ipair2 = set(ipair_to_pair2)

    ipair_only1 = ipair1 - ipair2
    ipair_only2 = ipair2 - ipair1

    if opts.fill_pairs:
        all_ipair   = ipair1 | ipair2
        # ipair_to_pair
        ipair_to_pair = {}
        ipair_to_pair.update(ipair_to_pair1)
        ipair_to_pair.update(ipair_to_pair2)

        with open(tc_outfp, 'wb') as out_fd:
            fmt = "{:>12} {:>12} {:>20}\n"
            sorted_ipair = sorted(list(all_ipair))
            diff_iter = gen_diff_all(
                    sorted_ipair, ipair_to_tc1, ipair_to_tc2, opts.dval)
            for idon, iacc, tc in diff_iter:
                don, acc = ipair_to_pair[(idon,iacc)]
                out_fd.write(fmt.format(don, acc, tc))

    else:
        ok_ipair    = ipair1 & ipair2

        with open(tc_outfp, 'wb') as out_fd:
            fmt = "{:>12} {:>12} {:>20}\n"
            sorted_ipair = sorted(list(ok_ipair))
            diff_iter = gen_diff(sorted_ipair, ipair_to_tc1, ipair_to_tc2)
            for idon, iacc, tc in diff_iter:
                don, acc = ipair_to_pair1[(idon,iacc)]
                out_fd.write(fmt.format(don, acc, tc))

    fmt = "    {:>12} {:>12}"
    print("** Below pairs exist only in {}.".format(tc_fp1))
    for idon, iacc in sorted(list(ipair_only1)):
        don, acc = ipair_to_pair1[(idon,iacc)]
        print(fmt.format(don, acc))
    print("")

    print("** Below pairs exist only in {}.".format(tc_fp2))
    for idon, iacc in sorted(list(ipair_only2)):
        don, acc = ipair_to_pair2[(idon,iacc)]
        print(fmt.format(don, acc))
    print("")

    print("** finished completely. **")

def parse_options():
    """Parse and get the command line options."""

    # make argument parser
    import argparse
    parser = argparse.ArgumentParser(
        description= (
            'calculate difference between two energy conductivity files.\n'+
            'Result(that is written into output_tcdiff_file) = tc_file2 - tc_file1'
            ) )

    # add argument definitions
    parser.add_argument(
            'tc_fps', nargs='+', #action='append',
            help='specify filenames. ')

    parser.add_argument(
            '-o', '--output-filename', dest='output_fp', required=True,
            help='specify a filename to write the difference data.')

    parser.add_argument(
            '-f', '--fill-deficit-pairs', dest='fill_pairs', default=False,
            required=False, action='store_true',
            help='fill the deficit pairs.')

    parser.add_argument(
            '-d', '--deficit-value', dest='dval',
            type=float, required=False, default=10**(-8),
            help='the value to be set for deficit pairs if set -f.')

    # make arguments
    options = parser.parse_args()
    return options



def load_tcdata(fp):
    if fp.endswith('.gz'):
        import gzip
        fd = gzip.open(fp, 'rb')
    else:
        fd = open(fp, 'rb')

    lines = ( line for line in fd
            if not line.startswith('#')
            if not line.isspace() )

    for line in lines:
        cols = line.split()
        if len(cols) == 3:
            don, acc, tc = cols[0], cols[1], float(cols[2])
            rest = None
        elif len(cols) >= 4:
            don, acc, tc, rest = cols[0], cols[1], float(cols[2]), cols[3]

        yield don, acc, tc, rest

    fd.close()

def get_id(name):
    return name.split('_')[0]

def gen_diff(target_ipair, dict1, dict2):
    """Get dict2 - dict1."""
    for ipair in target_ipair:
        tc1 = dict1[ipair]
        tc2 = dict2[ipair]
        idon, iacc = ipair
        yield idon, iacc, tc2-tc1

def gen_diff_all(all_ipair, dict1, dict2, val):
    """Get dict2 - dict1 for all."""
    for ipair in all_ipair:
        tc1 = dict1.get(ipair, val)
        tc2 = dict2.get(ipair, val)
        idon, iacc = ipair
        yield idon, iacc, tc2-tc1


if __name__ == '__main__':
    main()
