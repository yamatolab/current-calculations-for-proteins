"""Tool console scripts."""

from __future__ import print_function

# Standard module
import argparse
import curp.tool.sele as sele
from curp.console import exec_command


manual = {
          'sele': ('Output an irEC file with removed or kept only selected '
                   'residues from input irEC file.'),
          'helix': '',
          'neighbor': '',
          'noasa': '',
          'nocap': 'Remove residue pairs that contain the capped residues.',
          'nohelix': '',
          'noloop': '',
          'noneighbor': ('Remove residue pairs that indicate covalent '
                         'peptide bonds.'),
          'threshold': ''
          }


def arg_sele(parser=None):

    # Initialize
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['sele']
    parser.prog = 'sele'

    # Add commands as subparsers
    sp = parser.add_subparsers(dest='sele', help='Curp tool select commands')

    cmd_list = manual.keys()
    cmd_list.remove('sele')
    cmd_sp = {cmd: sp.add_parser(cmd, help=manual[cmd]) for cmd in cmd_list}

    arg_sel_helix(cmd_sp['helix'])
    arg_sel_neighbor(cmd_sp['neighbor'])
    arg_sel_noasa(cmd_sp['noasa'])
    arg_sel_nocap(cmd_sp['nocap'])
    arg_sel_nohelix(cmd_sp['nohelix'])
    arg_sel_noloop(cmd_sp['noloop'])
    arg_sel_noneighbor(cmd_sp['noneighbor'])
    arg_sel_thresh(cmd_sp['threshold'])


    if return_parser:
        return parser


def arg_sel_helix(parser=None):
    # Initialize
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['helix']
    parser.prog = 'sel_helix'
    parser.set_defaults(func=sele.helix)

    parser.add_argument('helix_fn', help='')
    parser.add_argument('arg_3', default=None)
    parser.add_argument('-f','--ec-fn', dest='ec_fn', default=None,
                        help=('[Optional] irEC data filename. If not '
                              'provided requires standard input.'))

    if return_parser:
        return parser

def arg_sel_neighbor(parser=None):
    # Initialize
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['neighbor']
    parser.prog = 'sel_neighbor'
    parser.set_defaults(func=sele.neighbor)

    parser.add_argument('excluded_list', default=['WAT'], nargs='+')
    parser.add_argument('-f','--ec-fn', dest='ec_fn', default=None,
                        help=('[Optional] irEC data filename. If not '
                              'provided requires standard input.'))

    if return_parser:
        return parser


def arg_sel_noasa(parser=None):
    # Initialize
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['noasa']
    parser.prog = 'sel_noasa'
    parser.set_defaults(func=sele.noasa)

    parser.add_argument('asa_fn')
    parser.add_argument('-f','--ec-fn', dest='ec_fn', default=None,
                        help=('[Optional] irEC data filename. If not '
                              'provided requires standard input.'))

    if return_parser:
        return parser


def arg_sel_nocap(parser=None):
    # Initialize
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['nocap']
    parser.prog = 'sel_nocap'
    parser.set_defaults(func=sele.nocap)

    parser.add_argument('-f','--ec-fn', dest='ec_fn', default=None,
                        help=('[Optional] irEC data filename. If not '
                              'provided requires standard input.'))

    if return_parser:
        return parser


def arg_sel_nohelix(parser=None):
    # Initialize
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['nohelix']
    parser.prog = 'sel_nohelix'
    parser.set_defaults(func=sele.nohelix)

    parser.add_argument('-f','--ec-fn', dest='ec_fn', default=None,
                        help=('[Optional] irEC data filename. If not '
                              'provided requires standard input.'))
    parser.add_argument('helix_fn')
    parser.add_argument('arg_3', default=None)

    if return_parser:
        return parser


def arg_sel_noloop(parser=None):
    # Initialize
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['noloop']
    parser.prog = 'sel_noloop'
    parser.set_defaults(func=sele.noloop)

    parser.add_argument('ss_fn')
    parser.add_argument('excluded_list', nargs='+', default=['WAT'])
    parser.add_argument('-f','--ec-fn', dest='ec_fn', default=None,
                        help=('[Optional] irEC data filename. If not '
                              'provided requires standard input.'))

    if return_parser:
        return parser


def arg_sel_noneighbor(parser=None):
    # Initialize
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['noneighbor']
    parser.prog = 'sel_noneighbor'
    parser.set_defaults(func=sele.noneighbor)

    parser.add_argument('excluded_list', nargs='+', default=['WAT'])
    parser.add_argument('-f','--ec-fn', dest='ec_fn', default=None,
                        help=('[Optional] irEC data filename. If not '
                              'provided requires standard input.'))

    if return_parser:
        return parser


def arg_sel_thresh(parser=None):
    # Initialize
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['threshold']
    parser.prog = 'sel_thresh'
    parser.set_defaults(func=sele.threshold)

    parser.add_argument('-f','--ec-fn', dest='ec_fn', default=None,
                        help=('[Optional] irEC data filename. If not '
                              'provided requires standard input.'))


def main():
    parser = arg_sele()
    exec_command(parser)


if __name__ == '__main__':
    main()
