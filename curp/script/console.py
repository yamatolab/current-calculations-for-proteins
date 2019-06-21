"""Script console arguments."""

from __future__ import print_function

#Standard module
import argparse

import curp


manual = {
          'conv-trj': 'Various scripts to process and analyze trajectories.',
          'analyze': 'Various scripts to analyze CURP results.',
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


