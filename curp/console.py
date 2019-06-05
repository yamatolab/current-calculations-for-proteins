"""Console scripts.

Argument parsing for the console use of curp."""

from __future__ import print_function

# Standard module
import argparse

import curp

manual = {
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

          'conv-trj': 'Various scripts to process and analyze trajectories.',

          'cal-tc': 'Calculate transport coefficients from energy flux data.',

          'graph-een': 'Show network chart of the energy conductivities.',

          'analyze': 'Various scripts to analyze CURP results.',

          'sum-acf': ('Average auto-correlation function over the given '
                      'trajectories'),

          # Descriptions used in analyze
          'divide-flux': ('Divide the flux data of all into the flux data for '
                          'each donor and acceptor.'),

          'pickup-respairs': ('Pick up the residue pair table to given file '
                              'as a ndx format.'),

          'simplify-tensor': 'Show and make figure for stress ratio.',

          'sum-tc': 'Summarize transport coefficient',

          'get-ncdata': ('Get simple text data from file in netcdf format '
                         'by given arguments.')
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
    sp_conv_trj = sp.add_parser('conv-trj',
                                help=manual['conv-trj'])
    sp_analyze = sp.add_parser('analyze',
                               help=manual['analyze'])

    # Add arguments and subparsers to each command.
    arg_compute(sp_compute)
    arg_cal_tc(sp_cal_tc)
    arg_graph_een(sp_graph_een)
    arg_sum_acf(sp_sum_acf)
    arg_conv_trj(sp_conv_trj)
    arg_analyze(sp_analyze)

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


def top_n_trj(func):
    """Decorator to make topology and trajectory objects before the function.

    Most functions called by conv-trj command need coordinate and trajectory
    objects. This decorator include their creation before calling the function.

    Parameters
    ----------
    func: function or method
    """

    def new_func(**kwds):

        if kwds['is_crd']:
            kwds['trj_type'] = 'crd'
        elif kwds['is_vel']:
            kwds['trj_type'] = 'vel'

        # Get topology and parameter file.
        tpl = curp.get_tpl(kwds['tpl_fmt'], kwds['tpl_fp'])
        # Get trajectory.
        trj = curp.gen_trj(tpl, **kwds)

        func(tpl, trj, **kwds)

    return new_func


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


def arg_conv_trj(parser=None):
    """Parse command line arguments of curp conv-trj"""

    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['conv-trj']
    parser.prog = 'conv-trj'

    # Add subparsers
    sp = parser.add_subparsers(dest='proc_name',
                               help='Processing command help')

    convert_only = sp.add_parser('convert-only',
                                 help=('Convert the trajectory into'
                                       'other format.'))
    dryrun = sp.add_parser('dry-run',
                           help='Do dry-run mode.')
    adjust_vel = sp.add_parser('adjust-vel',
                               help=('Adjust time of the velocity '
                                     'trajectory to t from t-dt/2.'))
    mask = sp.add_parser('mask',
                         help='Remove solvent from trajectory.')
    dist = sp.add_parser('dist',
                         help='Calculate inter-residue distances.')

    # Set process parcers arguments
    # parser.set_defaults(func=...) associate a function to the subparser
    # top_n_trj is a decorator to get the topology and the trajectory.

    convert_only.set_defaults(
        func=top_n_trj(curp.script.do_convert_only))

    dryrun.set_defaults(
        func=top_n_trj(dryrun))

    adjust_vel.set_defaults(
            func=top_n_trj(curp.script.adjust_vel))

    mask.set_defaults(
            func=top_n_trj(curp.script.do_mask))
    mask.add_argument(
            '-m', '--mask-filename', metavar='MASK_FILENAME',
            dest='mask_fn', required=True,
            help=('The mask file name.'))

    dist.set_defaults(
            func=top_n_trj(curp.script.do_distance))
    dist.add_argument('-m', '--method', metavar='DISTANCE_METHOD',
                      dest='dist_method', required=True,
                      choices=['cog', 'nearest', 'farthest'],
                      help=('The method used in calculating '
                            'the inter-residue distances.')
                      )
    dist.add_argument('-c', '--cutoff-length', metavar='CUTOFF_LENGTH',
                      dest='dist_cutoff', required=True,
                      default='5.0', type=float,
                      help=('The distance cutoff in angstrom.'))

    # Add common arguments
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
            '-crd',  # metavar='TRJ_TYPE',
            dest='is_crd', required=False, action='store_true',
            help=("Specify the format of the trajectory."
                  "This argument allows you specify multiple formats."))

    group.add_argument(
            '-vel',  # metavar='TRJ_TYPE',
            dest='is_vel', required=False, action='store_true',
            help=("Specify the format of the trajectory."
                  "This argument allows you specify multiple formats."))

    # add argument definitions
    parser.add_argument(
            '-i', '--input-filenames', metavar='TRJ_FILENAME',
            dest='input_trj_fns', required=True, nargs='+', action='append',
            help=('The trajectory file names.'))

    parser.add_argument(
            '-if', '--input-formats', metavar='FORMAT',
            dest='input_trj_fmts', required=True, action='append',
            help=("Specify the format of the trajectory."
                  "This argument allows you specify multiple formats."))

    parser.add_argument(
            '--irange', metavar=('FIRST', 'LAST', 'INTER'),
            dest='input_fst_lst_int', required=False,
            default=None, type=int, nargs=3,  action='append',
            help='The trajectory range to process over all trajectory file.')

    parser.add_argument('-o', '--output-filename', metavar='TRJ_FILENAME',
                        dest='output_trj_fn', required=False,
                        help=('The trajectory file name.'))

    parser.add_argument('-of', '--output-format', metavar='FORMAT',
                        dest='output_trj_fmt', required=False,
                        help=('The trajectory file format for output.'))

    parser.add_argument('--orange', metavar=('FIRST', 'LAST', 'INTER'),
                        dest='output_fst_lst_int', required=False,
                        default=[0, -1, 1], type=int, nargs=3,
                        help=('The trajectory range for output'))

    parser.add_argument('-pf', '--topology-format', metavar='FORMAT',
                        dest='tpl_fmt', required=True,
                        help='Specify the format of the topology file.')

    parser.add_argument('-p', '--topology-file', metavar='TPL_FILENAME',
                        dest='tpl_fp', required=True,
                        help='Specify the topology file.')

    parser.add_argument(
            '-nf', '--name-format', dest='dist_format',
            required=False, default='{rid:05}_{rname}',
            help='Specify the format for representing residue identify.')

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


def arg_analyze(parser=None):
    """Parse command line arguments of curp analyze"""

    # Initialize
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['analyze']
    parser.prog = 'analyze'

    # Definitions
    sp = parser.add_subparsers(dest='process',
                               help=manual['analyze'])
    sp_divide = sp.add_parser('divide-flux',
                              help=manual['divide-flux'])
    sp_respairs = sp.add_parser('pickup-respairs',
                                help=manual['pickup-respairs'])
    sp_simplify = sp.add_parser('simplify-tensor',
                                help=manual['simplify-tensor'])
    sp_sum_tc = sp.add_parser('sum-tc',
                              help=manual['sum-tc'])
    sp_ncdata = sp.add_parser('get-ncdata',
                              help=manual['get-ncdata'])

    arg_divide_flux(sp_divide)
    arg_respairs(sp_respairs)
    arg_simplify(sp_simplify)
    arg_sum_tc(sp_sum_tc)
    arg_ncdata(sp_ncdata)

    if return_parser:
        return parser


# curp analyze subparsers
def arg_divide_flux(parser=None):
    """Parse command line arguments of curp analyze divide-flux"""
    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.description = manual['divide-flux']
    parser.proc = 'divide-flux'
    parser.set_defaults(func=curp.script.analyze.divide_flux)

    # Argument definition for divide
    parser.add_argument(
            'flux_fns', nargs='+',
            help=('Specify filenames in order. '
                  'ex.) flux.dat0001, flux.dat0002, ...'))
    parser.add_argument(
            '-o', '--output-filename', dest='output_fn', required=True,
            help='Specify a filename to output.')
    parser.add_argument(
            '-d', '--donor', dest='donor_line', required=False, default='',
            help='Donor names to output to files.')
    parser.add_argument(
            '-a', '--acceptor', dest='acceptor_line', required=False,
            default='',
            help='Acceptor names to output to files.')
    parser.add_argument(
            '-t', '--dt', dest='dt', required=True, default=1.0,
            help='Time width[ps] per 1 step.')
    parser.add_argument(
            '-c', '--column', dest='column_line', required=False, default='',
            help='Column numbers to output.')

    if return_parser:
        return parser


def arg_respairs(parser=None):
    """Parse command line arguments of curp analyze pickup-respairs"""

    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.definition = manual['pickup-respairs']
    parser.proc = 'pickup-respairs'
    parser.set_defaults(func=curp.script.analyze.pickup_respairs)

    parser.add_argument(
            dest='trj_fns', nargs='+', action='append',
            help=('The trajectory file names.'))
    parser.add_argument('-if', '--input-format', metavar='FORMAT',
                        dest='input_trj_fmt', required=True, action='append',
                        help='Specify formats of trajectories.')
    parser.add_argument(
            '-p', '--input-prmtop-file', dest='prmtop_fn', required=True,
            help='Specify the prmtop file to be amber format.')
    parser.add_argument(
            '-pf', '--input-prmtop-format', dest='prmtop_fmt', required=True,
            help='Specify the format of the prmtop file.')
    parser.add_argument(
            '-i', '--interval', dest='interval', default=1, type=int,
            help=('Specify interval step to perform the calculation for '
                  'trajectory.'))
    parser.add_argument(
            '-m', '--cutoff-method', dest='cutoff_method', default='nearest',
            required=False, choices=['com', 'nearest', 'farthest'],
            help='Cutoff method; com, nearest and farthest for residue.')
    parser.add_argument(
            '-c', '--cutoff', dest='cutoff',
            required=False, default=5.0,
            help='Specify the cutoff to pick up.')
    parser.add_argument(
            '-t', '--trim-resnames', dest='trim_resnames', nargs='*',
            required=False, default=[],
            help='Residue names for trimming from the target residues.')
    parser.add_argument(
            '-f', '--format', dest='format',
            required=False, default='{rid:05}_{rname}',
            help='Specify the format for representing residue identify.')
    parser.add_argument(
            '-U', '--disble-union', dest='is_union',
            required=False, action='store_false',
            help=('Whether residue pair table is calculated by union '
                  'or intersection when collecting over trajectoryspecify'))
    parser.add_argument('-e', '--exclude-resids', metavar='FIRST:LAST',
                        dest='ext_resids', required=False,
                        default='', nargs='*',
                        help=('Residues that you want to exclude.'
                              'ex) 5:10 80:150'))
    parser.add_argument(
            '-b', '--enable-residues-both-cxcluded', dest='is_both',
            required=False, action='store_true',
            help=('Whether include residue pairs that exclude residues'))

    if return_parser:
        return parser


def arg_simplify(parser=None):
    """Parse command line arguments of curp analyze simplify-tensor"""

    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.definition = manual['simplify-tensor']
    parser.proc = 'simplify-tensor'
    parser.set_defaults(func=curp.script.analyze.simplify_tensor)

    # add argument definitions
    parser.add_argument(
            '-i', '--input-data', dest='filename', required=True,
            help='Specify input filename for the stress data.')
    parser.add_argument(
            '-l', '--labels', dest='labels',
            required=False, default='',
            help=('Specify labels of components to analyze.'
                  'ex.) -l "total,bond,angle,..."'))
    parser.add_argument(
            'fns', nargs='*',    # action='append',
            help=('Specify additive filenames. '
                  'ex.) label_data1, label_data2, ...'))
    parser.add_argument(
            '-s', '--every-snapshot', dest='snapshot',
            required=False, action='store_true', default=False,
            help='Specify flag to average the magnitude for every snapshot.')

    if return_parser:
        return parser


def arg_sum_tc(parser=None):
    """Parse command line arguments of curp analyze sum-tc"""

    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.definition = manual['sum-tc']
    parser.proc = 'sum-tc'
    parser.set_defaults(func=curp.script.analyze.summarize_tc)

    # add argument definitions
    parser.add_argument('tc_fps', nargs='+',
                        help='Data files')

    if return_parser:
        return parser


def arg_ncdata(parser=None):
    """Parse command line arguments of curp analyze get-ncdata"""

    return_parser = parser is None
    if return_parser:
        parser = argparse.ArgumentParser()

    parser.definition = manual['get-ncdata']
    parser.proc = 'get-ncdata'
    parser.set_defaults(func=curp.script.analyze.get_ncdata)

    parser.add_argument(dest='acf_fp', metavar='ACF_FILE',
                        help='The filepath of auto-correlation function data.')
    parser.add_argument('-r', '--group-ranges', metavar='FRIST:LAST',
                        nargs='*',
                        dest='group_ranges', required=True,
                        default='',
                        help='The pair range list to want to gain.')
    parser.add_argument('-n', '--dataname',
                        dest='dataname', required=False,
                        default='acf',
                        help='The name of the netcdf data you want to gain.')
    parser.add_argument('-o', '--output-prefix',
                        dest='prefix', required=False,
                        default='acf',
                        help=('The prefix of the files to want to write, '
                              'that includes directory path.'))

    if return_parser:
        return parser


def main():
    parser = arg_curp()
    exec_command(parser)


if __name__ == '__main__':
    main()
