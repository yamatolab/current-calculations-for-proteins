import os, sys
curp_path = os.environ["CURP_HOME"]
src_path = os.path.join(curp_path, 'src')
sys.path.insert(0, src_path)

import parser
get_topology = parser.get_tplprm_simple

def get_topology_old(tplprm_fn):
    from amber.topology import TopologyParser, Format2AmberBaseConverter
    raw_tpl = TopologyParser(tplprm_fn)
    raw_tpl.parse()
    tpl = Format2AmberBaseConverter(raw_tpl, use_atomtype=False)
    tpl.convert()
    return tpl

def gen_crds(*args, **kwds): return gen_trj(trj_type='crd', *args, **kwds)

def gen_vels(*args, **kwds): return gen_trj(trj_type='vel', *args, **kwds)


def gen_trj(trj_fns, trj_fmt, trj_type=None,
        natom=None, fst_lst_int=(0,-1,1), use_pbc=True):

    if trj_type=='crd':
        Parser = parser.get_coordinate_parser(trj_fmt)
    elif trj_type=='vel':
        Parser = parser.get_velocity_parser(trj_fmt)
        use_pbc = False
    else:
        raise Exception()

    trj_pbcs = parser.gen_trajectory(trj_fns, Parser, natom, logger=None,
            use_pbc=use_pbc, first_last_interval=fst_lst_int)

    return trj_pbcs # yield istp, snap, box

def gen_trajectory_multi(trj_fns_list, trj_fmts, natom, fst_lst_int_list,
        trj_type='crd',logger=None, use_pbc=False):

    from itertools import izip_longest

    for trj_fns, fmt, fst_lst_int in izip_longest(
            trj_fns_list, trj_fmts, fst_lst_int_list):

        # determine parser
        Parser = get_any_parser(fmt, natom, use_pbc, trj_type)

        # generate trajectory 
        trj = parser.gen_trajectory(trj_fns, Parser, natom, logger=None, 
                use_pbc=use_pbc, first_last_interval=fst_lst_int)

        # iterate
        for istp, snap, box in trj:
            yield istp, snap, box

def get_any_parser(trj_fmt, natom, use_pbc, trj_type='crd'):

    if trj_fmt == 'restart':
        Parser = parser.get_parser(trj_fmt, 'restart')
        Parser.set_trjtype(trj_type)
        return Parser

    else:
        if trj_type == 'crd':
            return parser.get_parser(trj_fmt, 'coordinate')

        elif trj_type == 'vel':
            return parser.get_parser(trj_fmt, 'velocity')


def gen_trajectory_old(traj_fns, natom, interval=1, use_pbc=True):
    amber_parser_path = os.path.join(curp_path, 'src' , 'parser')
    sys.path.insert(0, amber_parser_path)
    from amber.trajectory import CoordinateParser
    import parser 
    crd_parser = parser.gen_trajectory(
            traj_fns, CoordinateParser, natom, use_pbc=use_pbc)

    for ntraj_1, crd, box in crd_parser:
        if (ntraj_1+1)%interval != 0: continue
        yield ntraj_1+1, crd, box

# table
from table.ini_parser import IniParser
from table.target import parse_target_atoms_line as parse_target_groups

# writer
from parser.writer import Writer as TrjWriter

import parallel
if parallel.use_mpi:
    ParallelProcessor = parallel.ParallelProcessor
else:
    ParallelProcessor = parallel.SequentialProcessor

