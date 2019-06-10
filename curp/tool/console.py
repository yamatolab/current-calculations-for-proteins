"""Tool console scripts."""

from __future__ import print_function

# Standard module
import argparse
import curp.tool as tool
from curp.console import exec_command


manual = {
          'tool': ''
          }


def arg_tool(parser=None):

    # Initialize
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['tool']
    parser.prog = 'tool'

    # Add commands as subparsers
    sp = parser.add_subparsers(dest='command', help='Curp tool commands')


def main():
    parser = arg_tool()
    exec_command(parser)


if __name__ == '__main__':
    main()
