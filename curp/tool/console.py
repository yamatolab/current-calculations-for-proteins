"""Tool console scripts."""

from __future__ import print_function

# Standard module
import argparse
import curp.tool as tool
from curp.console import exec_command


manual = {
          'tool': 'Scripts to process irEC data',
          'renum-residue': 'Renumbering the residue numbers'
          }


def arg_tool(parser=None):
    # Initialize
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['tool']
    parser.prog = 'tool'

    # Add commands from manual dictionnary as subparsers
    sp = parser.add_subparsers(dest='command', help='Curp tool commands')

    cmd_list = manual.keys()
    cmd_list.remove('tool')
    cmd_sp = {cmd: sp.add_parser(cmd, help=manual[cmd]) for cmd in cmd_list}

    # Add deeper subparsers
    from curp.tool.sele.console import manual as sele_manual
    cmd_sp['sele'] = sp.add_parser('sele',
                                   help=sele_manual['sele'])

    # Configure subparsers
    arg_renum_res(cmd_sp['renum-residue'])

    from curp.tool.sele.console import arg_sele
    arg_sele(cmd_sp['sele'])


def arg_renum_res(parser=None):
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['renum-residue']
    parser.prog = 'renum-residue'
    parser.set_defaults(func=tool.renum_residue.main)

    # Check if there is no input
    # from sys import stdin
    # if stdin.isatty():
    #     print('no stdin')
    #     parser.add_argument('ec_fn',
    #                         help=('Interresidue energy conductivity '
    #                               '(irEC) filename.'))
    # else:
    #     parser.set_defaults(ec_file=stdin)

    parser.add_argument('renums', nargs='+',
                        help=('Residue to renumber and which number to assign.'
                              'Ex: 2:10 would add 9 to every residue numbers, '
                              'starting at residue number 2.'))
    parser.add_argument('-f', '--ec-fn', dest='ec_fn', default=None,
                        help=('[Optional] irEC data filename. If not '
                              'provided requires standard input.'))

    if return_parser:
        return parser


def main():
    parser = arg_tool()
    exec_command(parser)


if __name__ == '__main__':
    main()
