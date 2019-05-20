# Getting main curp function.
from .curp import curp

# Topology
from.parser import get_tplprm_simple as get_tpl

# Trajectory
from .script.conv_trj import gen_trj

# Get all the nice scripts
import .script as script


# writer
from .parser.writer import Writer as TrjWriter

import curp.parallel as parallel
if parallel.use_mpi:
    ParallelProcessor = parallel.ParallelProcessor
else:
    ParallelProcessor = parallel.SequentialProcessor
