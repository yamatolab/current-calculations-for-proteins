"""Console scripts.

Argument parsing for the console use of curp."""

from __future__ import print_function

# Standard module
import argparse

import curp

manual = {
          'cal-tc': 'Calculate transport coefficients from energy flux data.',

          'curp': ("CURP - CURrent calculation for Proteins.\n"
                   "Compute inter-residue flow of energy or heat and atomic"
                   "stress tensors in a protein. Create Inter-residue Energy"
                   "Exchange Network as a graph.\n\n"
                   "This package also contains functionalities to format "
                   "the trajectories and analyze said results.\n\n"
                   "Developed by Yamato's Laboratory, University of Nagoya, "
                   "Japan."),

          'compute': ('Launch flux or stress tensor computations given a '
                      'configuration file.'),

          'graph-een': 'Show network chart of the energy conductivities.',

          'sum-acf': ('Average auto-correlation function over the given '
                      'trajectories'),
        }


def arg_curp():
    """Parse command line arguments of curp command."""

    parser = argparse.ArgumentParser(description=manual['curp'],
                                     prog='curp')

    # Add commands as subparsers.
    sp = parser.add_subparsers(dest='command', help='Curp commands')
    sp_compute = sp.add_parser('compute',
                               help=manual['compute'])
    sp_cal_tc = sp.add_parser('cal-tc',
                              help=manual['cal-tc'])
    sp_graph_een = sp.add_parser('graph-een',
                                 help=manual['graph-een'])
    sp_sum_acf = sp.add_parser('sum-acf',
                               help=manual['sum-acf'])

    from curp.script.console import manual as script_manual
    sp_conv_trj = sp.add_parser('conv-trj',
                                help=script_manual['conv-trj'])
    sp_analyze = sp.add_parser('analyze',
                               help=script_manual['analyze'])

    from curp.tool.console import manual as tool_manual
    sp_tool = sp.add_parser('tool',
                            help=tool_manual['tool'])

    # Add arguments and subparsers to each command.
    arg_compute(sp_compute)
    arg_cal_tc(sp_cal_tc)
    arg_graph_een(sp_graph_een)
    arg_sum_acf(sp_sum_acf)

    from curp.script.console import arg_conv_trj, arg_analyze
    arg_conv_trj(sp_conv_trj)
    arg_analyze(sp_analyze)

    from curp.tool.console import arg_tool
    arg_tool(sp_tool)

    return(parser)


def exec_command(parser):
    """Call scripts depending on command line argument.

    Parameters
    ----------
    parser : argparse.ArgumentParser
    """

    args = parser.parse_args()
    command = args.func
    command(**vars(args))


def arg_compute(parser=None):
    """Parse command line arguments of curp compute"""

    # Initialize
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['compute']
    parser.prog = 'compute'

    # Set the default function used by the parser.
    parser.set_defaults(func=curp.curp)

    # Definitions
    parser.add_argument('-v', '--vervose',
                        dest='vervose', default=False,
                        action='store_true',
                        help='Print out informations.')

    parser.add_argument('-s', '--enable-serial',
                        dest='use_serial', default=False,
                        action='store_true',
                        help=("Calculate in serial, don't calculate"
                              "in parallel."))

    parser.add_argument('--output-conf-default',
                        dest='output_conf_def', default=False,
                        action='store_true',
                        help=('Output the config parameters in ini format '
                              'with default values.'))

    parser.add_argument('--output-conf-formatted',
                        dest='output_conf_fmtd', default=False,
                        action='store_true',
                        help=('Output the config parameters in rest style '
                              'with default values.'))

    parser.add_argument('input_', nargs='?',
                        default='run.cfg',
                        help='Specify input filenames.')

    if return_parser:
        return parser


def arg_cal_tc(parser=None):
    """Parse the command line arguments of curp cal-tc."""

    # Make argument parser
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['cal-tc']
    parser.prog = 'cal-tc'
    parser.set_defaults(func=curp.script.cal_tc)

    # add argument definitions
    parser.add_argument(dest='flux_fn', metavar='FLUX_FILENAME',
                        help=('The filepath of flux data.'))

    parser.add_argument('-o', '--tc-file', metavar='FILE',
                        dest='tc_fn', required=True,
                        help=('The filename of tc data.'))

    parser.add_argument('-a', '--acf-file', metavar='FILE',
                        dest='acf_fn', required=False,
                        default='',
                        help=('The filename of acf data.'))

    parser.add_argument('-af', '--acf-file-format', metavar='FORMAT',
                        dest='acf_fmt', required=False,
                        default="netcdf", choices=['ascii', 'netcdf'],
                        help=('The file format of acf data.'))

    parser.add_argument('-r', '--frame-range',
                        metavar=('FIRST', 'LAST', 'INTER'),
                        dest='frame_range', required=False,
                        default=[1, -1, 1], type=int, nargs=3,
                        help='The range of frame; first frame, last frame, '
                              'interval step.')

    parser.add_argument('-s', '--average-shift', metavar='SHIFT',
                        dest='avg_shift', required=False,
                        default=1, type=int,
                        help='The frame to shift for averaging.')

    parser.add_argument('--sample-number', dest='nsample',
                        type=int, required=False, default=0,
                        help=('number of sample for one flux data file. '
                              'Default value is 0 that present to make samples'
                              ' as much as possible'))

    parser.add_argument(
        '-dt', '--dt', metavar='DT', dest='d_t', required=False,
        default=None, type=float,
        help=('t of between the neighbour frames. The unit is in ps. '
              + 'Default value is determined by time variable in flux data.'))

    parser.add_argument('-c', '--coefficient', metavar='COEFFICIENT',
                        dest='coef', required=False,
                        default=1.0, type=float,
                        help='Multiply acf by given coefficient.')

    parser.add_argument('-v', '--vervose', dest='use_debug',
                        action='store_true', required=False,
                        help='turn on debug mode.')

    parser.add_argument('--no-axes', dest='no_axes',
                        action='store_true', required=False,
                        help='with this option, scalar flux is handled.')

    if return_parser:
        return parser


def arg_graph_een(parser=None):
    """Parse the command line arguments of curp graph-een."""

    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['graph-een']
    parser.proc = 'graph-een'
    parser.set_defaults(func=curp.script.graph_een)

    # Add argument definitions
    parser.add_argument(
            'data_fns', nargs='+',
            help='Specify filenames. ')
    parser.add_argument(
            '-t', '--targets', dest='targets', required=False,
            nargs='*', default=['1-'],
            help=('This option allow you to show only target '
                  'and thier neighbhor nodes. '
                  'ex), [-t 1-5], [-t 2-9 15 17], [-t 25-].'))
    parser.add_argument(
            '--forced-output-nodes', dest='force_nodes', required=False,
            nargs='*', default=[],
            help=('nodes to forcibly show.'
                  'ex), [-t 1-5], [-t 2-9 15 17], [-t 25-].'))
    parser.add_argument(
            '-p', '--bring-node-pair-close-together', dest='close_pairs',
            required=False, nargs='*', default=[],
            help=('bring the node pair close together on the EEN graph.'
                  'ex), 1:2, 3:5, 3:15.'))
    parser.add_argument(
            '-c', '--cluster-filename', dest='cluster_fn',
            required=False, default='',
            help=('Specify a cluster file name. '))
    parser.add_argument('-n', '--node-style-filename', dest='node_fn',
                        required=False, default='',
                        help=('Specify a file name which node '
                              'style difinitions are written.')
                        )
    parser.add_argument(
            '--with-one-letter', dest='use_1letter', default=False,
            required=False, action='store_true',
            help='use 1 letter representation for amino acid name.')
    parser.add_argument(
            '-f', '--output-een-filename',
            dest='fig_fn', required=False,
            help='Specify a filename to graph EEN.')
    parser.add_argument(
            '-r', '--threshold', dest='threshold',
            type=float, required=False, default=None,
            help='threshold value to show nodes on chart figure.')
    parser.add_argument(
            '-R', '--topology-threshold', dest='tpl_threshold',
            type=float, required=False, default=None,
            help=("threshold value to define graph topology. "
                  "If you don't give this parameter, "
                  "the value of threshold to show nodes will be used."))

    # Esthetic parameters
    parser.add_argument(
            '--ratio', dest='ratio',
            type=float, required=False, default=None,
            help=("ratio of height/width to chart a figure. "
                  "If '--graph-size is set, this option will be overwritten"
                  "by it."))
    parser.add_argument(
            "-s", "--graph-size", dest="graph_size",
            required=False, default=None,
            help=("The graph size in inch unit. ex) -s '3.4,2.5'"))
    parser.add_argument(
            '--title', dest='title', required=False, default='',
            help='Specify the title of the figure.')
    parser.add_argument(
            '--direction', dest='direction',
            required=False, choices=['LR', 'TB'], default='TB',
            help='Specify the kind of the direction.')
    parser.add_argument(
            '-I', '--hide-isolated-nodes', dest='hide_isonode',
            required=False, action='store_true',
            help='Hide isolated nodes when appling multiple irEC files.')
    parser.add_argument(
            '--show-negative-values', dest='use_decrease',
            action='store_true',
            help="show the nodes pair which have negative values.")
    parser.add_argument(
            '--alpha', dest='alpha',
            type=int, required=False, default=None,
            help='The transparency value of the background color.')
    parser.add_argument(
            '-lv', '--line-values', dest='line_values', required=False,
            nargs='*', type=float, default=[0.015, 0.008, 0.003],
            help=('The threshold values for line attributes. '
                  'Number of elements in the list must be equal '
                  'with all line attribute.'))
    parser.add_argument(
            '-lc', '--line-colors', dest='line_colors', required=False,
            nargs='*', type=str, default=['red', 'blue', 'green'],
            help=('The colors of line. Number of elements in the list '
                  'must be equal with all line attribute.'))
    parser.add_argument(
            '-lt', '--line-thicks', dest='line_thicks', required=False,
            nargs='*', type=float, default=[4.0, 4.0, 2.5],
            help=('The thickness of line. Number of elements in the list '
                  'must be equal with all line attribute.'))
    parser.add_argument(
            '-lw', '--line-weights', dest='line_weights', required=False,
            nargs='*', type=float, default=[5.0, 3.0, 1.0],
            help=('The weight of line. Number of elements in the list must'
                  'be equal with all line attribute.'))

    if return_parser:
        return parser


def arg_sum_acf(parser=None):
    """Parse command line arguments of curp sum-acf"""

    # Initialize
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['sum-acf']
    parser.prog = 'sum-acf'

    # Set the default function used by the parser.
    parser.set_defaults(func=curp.script.sum_acf)

    # Definitions
    parser.add_argument(
            dest='acf_fps', metavar='ACF_FILE', nargs='+',
            help=('The filepath of auto-correlation function data.'))
    parser.add_argument('-a', '--acf-file', metavar='FILE',
                        dest='output_acf_fp', required=True,
                        default='',
                        help=('The filepath of acf data.'))
    parser.add_argument('-t', '--tcs-file', metavar='FILE',
                        dest='tcs_fp', required=False,
                        default='',
                        help=('The filepath of tc time series data.'))
    parser.add_argument('-c', '--coefficient', metavar='COEFFICIENT',
                        dest='coef', required=False,
                        default=1.0, type=float,
                        help='Multiply acf by given coefficient.')

    if return_parser:
        return parser


def main():
    parser = arg_curp()
    exec_command(parser)


if __name__ == '__main__':
    main()
