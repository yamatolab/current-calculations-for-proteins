from __future__ import print_function

# standard modules
import os, sys

# curp modules
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path: sys.path.insert(0, topdir)
from exception import CurpException

# exceptions
class TrajectoryExistsError(CurpException): pass


class Writer:

    def __init__(self, dynamics_setting):

        dyn_sets = dynamics_setting

        crd_freq = dyn_sets.crds_frequency
        vel_freq = dyn_sets.vels_frequency
        last = 10**18

        crds_range = (crd_freq, last, crd_freq)
        vels_range = (crd_freq, last, vel_freq)

        import parser.writer as pw

        # coordinate
        if dyn_sets.crds_frequency >= 1:
            self.__crds_writer = pw.Writer(
                    trj_fp      = dyn_sets.crds_file[0] , 
                    trj_fmt     = dyn_sets.trj_format   , 
                    dt          = dyn_sets.dt           , 
                    is_vel      = False                 , 
                    fst_lst_int = crds_range)
        else:
            self.__crds_writer = None

        # velocity
        if dyn_sets.crds_frequency >= 1:
            self.__vels_writer = pw.Writer(
                    trj_fp      = dyn_sets.vels_file[0] , 
                    trj_fmt     = dyn_sets.trj_format   , 
                    dt          = dyn_sets.dt           , 
                    is_vel      = True                  , 
                    fst_lst_int = vels_range)

        else:
            self.__vels_writer = None

    def write_header(self):

        if self.__crds_writer:
            self.__crds_writer.write_header()

        if self.__vels_writer:
            self.__vels_writer.write_header()

    def write(self, istp, crd=None, vel=None):
        """Write trajectory."""

        if self.__crds_writer:
            self.__crds_writer.write(istp, crd)

            # The box is not implemented in dynamics run.
            #self.__crds_writer.write(istp, crd, box)

        if self.__vels_writer:
            self.__vels_writer.write(istp, vel)

