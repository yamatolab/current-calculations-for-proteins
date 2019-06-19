#! /usr/bin/env python2
from __future__ import print_function

import sys


def parse_line(line):
    cols = line.split()
    don_line, acc_line, rest = cols[0], cols[1], cols[2:]

    rid_don, rname_don = don_line.split('_')
    rid_don = int(rid_don)
    rid_acc, rname_acc = acc_line.split('_')
    rid_acc = int(rid_acc)

    return rid_don, rid_acc, rname_don, rname_acc, rest


def main(renums, ec_fn=None, **kwds):
    """Renumber the residue numbers according to renums.

    Parameters
    ----------
    renums : iterable of strings
        Residue to renumber and which number to assign.

    Examples
    --------
    With an irEC file ec.dat so that:
    >>> with open('ec.dat', 'r') as ec_file:
    ...     print(ec_file.read())
    ...
        00001_ACE    00002_ILE  0.1  0.1
        00001_ACE    00004_PRO  0.2  0.2
        00004_PRO    00006_GLU  0.3  0.3

    Then to renumber the first residue to 3, and the fourth to 10:
    >>> main('ec.dat', ['1:3', '4:10'])
        00003_ACE    00004_ILE  0.1  0.1
        00003_ACE    00010_PRO  0.2  0.2
        00010_PRO    00012_GLU  0.3  0.3
    """

    # Remove the only space line and the comment line

    if ec_fn is None:
        ec_file = sys.stdin
    else:
        ec_file = open(ec_fn, 'r')
    lines = [line.strip() for line in ec_file
             if not line.isspace() if not line.startswith('#')]

    # Obtain maximum residue number
    rid_max = 0
    for line in lines:
        rid1, rid2, rname1, rname2, rest = parse_line(line)
        rid_max = max(rid_max, rid1, rid2)

    new_rids = range(rid_max+1)
    old_rids = range(rid_max+1)

    # Make dictionary to alter the given residue identifiers
    for renum in renums:
        rid_old, rid_new = renum.split(':')
        rid_old, rid_new = int(rid_old), int(rid_new)

        shift = rid_new - rid_old
        new_rids[rid_old:] = [rid+shift for rid in old_rids[rid_old:]]

    # Convert residue identifiers into new ones
    for line in lines:
        rid1, rid2, rname1, rname2, rest = parse_line(line)
        ec = rest[0]

        fmt = '{:05}_{:<7}  '
        print(fmt.format(new_rids[rid1], rname1),
              fmt.format(new_rids[rid2], rname2), ec)


if __name__ == '__main__':
    from curp.tool.console import exec_command, arg_renum_residue
    parser = arg_renum_residue()
    exec_command(parser)
