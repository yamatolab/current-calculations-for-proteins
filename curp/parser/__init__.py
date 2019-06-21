"""
Gives functions to create topology and trajectory variables from files indicated in .cfg configuration file.
"""

import os, sys
import traceback
import itertools as it
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)

################################################################################
from exception import CurpException
class FormatNotFoundError(CurpException): pass
class ParserNotFoundError(CurpException): pass
class ConverterNotFoundError(CurpException): pass
class NotFoundDictionary(CurpException): pass
class StepNumberInvalidError(CurpException): pass

import trans_rot

################################################################################
def get_tplprm(file_format, format_setting, potential_type, use_atomtype):

    # decide topology parser
    TopologyParser = get_topology_parser(file_format)
    tpl_parser = TopologyParser( format_setting.topology_file[0] )
    tpl_parser.parse()

    # convert topology to forcefield
    TopologyConverter = get_topology_converter(file_format, potential_type)
    tpl_conv = TopologyConverter(tpl_parser, use_atomtype)
    tpl_conv.convert(format_setting)

    return tpl_conv

def get_tplprm_simple(file_format, fp, use_atomtype=False):

    # decide topology parser
    TopologyParser = get_topology_parser(file_format)
    tpl_parser = TopologyParser( fp )
    tpl_parser.parse()

    # convert topology to forcefield
    TopologyConverter = get_topology_converter(file_format, file_format)
    tpl_conv = TopologyConverter(tpl_parser, use_atomtype)
    tpl_conv.convert()

    return tpl_conv


def gen_phase_trajectory(file_format, fmt_section,
                         first_last_interval, natom, use_pbc,
                         rem_trans, rem_rotate, target_atoms, masses,
                         logger):

    # define the parser class
    if fmt_section.target == 'restart':

        rst_format = fmt_section.restart_format
        rst_filname = fmt_section.restart_file[0]

        RestartParser = get_restart_parser(rst_format)

        rst_parser = RestartParser(rst_filname, natom)
        crd, vel, pdc = rst_parser.parse()

        # remove translation and rotation
        new_crd, new_vel = trans_rot.get_crd_vel_trans_rot_removed(
                masses, crd, vel, target_atoms, rem_trans, rem_rotate)
        crd_pbcs = [(0, new_crd, None)]
        vel_pbcs = [(0, new_vel, None)]

    elif fmt_section.target == 'trajectory':

        crd_format = fmt_section.coordinate_format
        vel_format = fmt_section.velocity_format
        crd_filename = fmt_section.coordinate_file
        vel_filename = fmt_section.velocity_file

        CoordinateParser = get_coordinate_parser(crd_format)
        VelocityParser   = get_velocity_parser(vel_format)

        crd_pbcs = gen_trajectory(fmt_section.coordinate_file,
                CoordinateParser, natom, logger, use_pbc, 
                first_last_interval=first_last_interval)

        vel_pbcs = gen_trajectory(fmt_section.velocity_file,
                VelocityParser, natom, logger, False,
                first_last_interval=first_last_interval)

    else:
        pass

    first, last, inter = first_last_interval
    logger.info_title('The range of trajectories')
    logger.info('first, last, interval of step = {} {} {}'.format(
        first, last, inter))
    yield

    istep = 0
    for (crd_data, vel_data) in it.izip(crd_pbcs, vel_pbcs):
        cstep_crd, crd, pbc_crd = crd_data
        cstep_vel, vel, pbc_vel = vel_data
        istep += 1
        logger.set_curstep(istep)

        if cstep_crd != cstep_vel:
            msg = 'crd step = {}, but vel step = {}.'
            raise StepNumberInvalidError( msg.format(cstep_crd, cstep_vel) )

        logger.debug_cycle('\n\n' + 20*'=' + 'STEP {}'.format(istep) + 20*'=')
        logger.debug_cycle('    current step = {}'.format(cstep_crd))
        logger.info_cycle('*** ISTEP = {}, CURRENT STEP = {} ***'
                .format(istep, cstep_crd))
        logger.debug_cycle()

        # remove translation and rotation
        new_crd, new_vel = trans_rot.get_crd_vel_trans_rot_removed(
                masses, crd, vel, target_atoms, rem_trans, rem_rotate)

        yield cstep_crd, (new_crd, new_vel, pbc_crd)

def gen_matrix_ensemble(file_format, fmt_section, first_last_interval, natom):
    # define the parser class

    return istep, matrix

################################################################################

def get_topology_parser(format):
    return get_parser(format, 'topology')

def get_coordinate_parser(format):
    return get_parser(format, 'coordinate')

def get_velocity_parser(format):
    return get_parser(format, 'velocity')

def get_restart_parser(format):
    return get_parser(format, 'restart')

def get_topology_converter(format, forcefield):
    dictname = 'converter_dict'
    key = (format, forcefield)

    # get the parser's directory
    parser_path = os.path.dirname(__file__)

    # get the plugin names
    modpaths = []
    for fd in os.listdir(parser_path):
        abspath = os.path.join(parser_path, fd)
        if os.path.isdir(abspath) and not fd.startswith('_'):
            modpaths.append(abspath)

    for mpath in modpaths:
        modname = os.path.split(mpath)[-1]
        try:
            module = load_module(modname, parser_path)
            if hasattr(module, dictname):
                dictionary = getattr(module, dictname)
            else:
                raise NotFoundDictionary(dictname)

            if key in dictionary:
                Converter = dictionary[key]
                break

        except ImportError:
            print(traceback.format_exc())
            # raise FormatNotFoundError(format)

        except:
            raise

    else:
        raise ParserNotFoundError(key)

    return Converter


def gen_trajectory(traj_files, Parser, natom, logger=None,
        use_pbc=False, first_last_interval=(1, -1, 1)):
    """Generate trajectory from files."""

    # logger
    class Logger:
        def info(self, *args, **kwds):
            pass

    if logger is None:
        logger = Logger()

    # range_line
    first, last, inter = first_last_interval

    int_count = 0
    istep_all = 0
    last = 10**18 if last <= -1 else last

    # do loop the trajectory files
    for traj_fn in traj_files:
        logger.info('*** reading trajectory file ***')
        logger.info('    {}'.format(traj_fn) )
        logger.info()
        parser = Parser(traj_fn, natom, use_pbc=use_pbc)
        # do loop the trajectory in a file
        for istep_1, traj in enumerate(parser):
            int_count += 1
            istep = istep_1 + 1
            istep_all += 1
            body, pbc = traj

            if istep_all < first:
                pass
            elif istep_all == first:
                yield istep_all, body, pbc
                int_count = 0
            elif (first < istep_all <= last) and (int_count == inter):
                int_count = 0
                yield istep_all, body, pbc
            elif last < istep_all:
                break
            else:
                continue

        parser.close()

def gen_matrix(matrix_files, Parser, natom, logger):
    pass


################################################################################
def get_parser(format, parser_name):
    dictname = parser_name + '_dict'

    # get the parser's directory
    parser_path = os.path.dirname(__file__)

    # get the plugin names
    modpaths = []
    for fd in os.listdir(parser_path):
        abspath = os.path.join(parser_path, fd)
        if os.path.isdir(abspath) and not fd.startswith('_'):
            modpaths.append(abspath)

    for mpath in modpaths:
        modname = os.path.split(mpath)[-1]
        try:

            # get parser dictionary
            module = load_module(modname, parser_path)
            if hasattr(module, dictname):
                dictionary = getattr(module, dictname)
            else:
                raise NotFoundDictionary(dictname)

            # get Parser
            if format in dictionary:
                Parser = dictionary[format]
                break

        except ImportError:
            print(traceback.format_exc())
            # raise FormatNotFoundError(format)

        except:
            raise

    else:
        raise ParserNotFoundError(format)

    return Parser


import imp
def load_module(modname,  basepath):
    f,n,d = imp.find_module(modname, [basepath])
    return imp.load_module(modname, f, n, d)


if __name__ == '__main__':
    TopologyParser = get_topology_parser('presto')
    tparser = TopologyParser('./presto_v2/test/ala_all.tpl')
    tparser.parse()

    Converter = get_topology_converter('presto', 'amber99')
    con = Converter(tparser)
    con.convert()
