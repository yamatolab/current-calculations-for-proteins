from __future__ import print_function

import os, sys
import time
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)
import utility
import numpy

# curp modules
import clog as logger

import base

################################################################################
class StressCurrentCalculator(base.CurrentCalculator):

    """
    A class to caluculate the stress current.
    """

    def __init__(self):
        base.CurrentCalculator.__init__(self)

    def prepare(self, *args, **kwds):
        """Prepare the calculations of all snapshots."""
        base.CurrentCalculator.prepare(self, *args, **kwds)
        
        self.scal = StressTensor( self.get_target_atoms(),
                self.get_iatm_to_igrp(), self.get_bonded_pairs())

    def cal_kinetic(self, vel):

        masses   = self.get_topology().get_atom_info()['masses']
        volumes  = self.get_volume()
        gvolumes = self.get_gvolume()

        t0 = time.time()

        current, current_grp, current_grp_inner = self.scal.cal_kinetic(
                vel, masses, volumes, gvolumes)

        t1 = time.time()
        self.store_time('kinetic current', t1-t0)

        return current, current_grp, current_grp_inner

    def cal_bonded(self, crd, bond_type):
        """Calculate the current for one of bonded types."""
        tbcal = self.get_tbforce()

        # calculate two-body force of each bonded interactions
        t0 = time.time()
        tbfs_info = tbcal.cal_bonded(bond_type)
        tbfs      = tbfs_info['tbforces']

        # get valume
        volumes   = self.get_volume()
        gvolumes  = self.get_gvolume()

        # calculate current
        t1 = time.time()
        current, current_grp, current_grp_inner = self.scal.cal_bonded(
                crd, tbfs, volumes, gvolumes)

        # log the elasped times
        t2 = time.time()
        self.store_time(bond_type+' pairwise', t1-t0)
        self.store_time(bond_type+' current' , t2-t1)

        return current, current_grp, current_grp_inner

    def cal_coulomb(self, crd):
        """Calculate the current for coulomb components."""

        table     = self.get_interact_table()
        type_func = self.get_tbforce().cal_coulomb
        cutoff    = self.get_setting().curp.coulomb_cutoff_length
        volumes   = self.get_volume()
        gvolumes  = self.get_gvolume()

        t0 = time.time()
        gen_tbfs = ( type_func(t)['tbforces'] for t in table )
        current, current_grp, current_grp_inner  = self.scal.cal_nonbonded(
                crd, gen_tbfs, volumes, table, gvolumes)

        t1 = time.time()
        dt_current = self.scal.dt

        self.store_time('coulomb pairwise' , t1 - t0 - dt_current)
        self.store_time('coulomb current' , dt_current)

        return current, current_grp, current_grp_inner

    def cal_vdw(self, crd):
        """Calculate the current for vdW components."""

        table     = self.get_interact_table()
        type_func = self.get_tbforce().cal_vdw
        cutoff    = self.get_setting().curp.vdw_cutoff_length
        volumes   = self.get_volume()
        gvolumes  = self.get_gvolume()

        t0 = time.time()

        gen_tbfs = ( type_func(t)['tbforces'] for t in table )
        current, current_grp, current_grp_inner  = self.scal.cal_nonbonded(
                crd, gen_tbfs, volumes, table, gvolumes)

        t1 = time.time()
        dt_current = self.scal.dt

        self.store_time('vdw pairwise' , t1 - t0 - dt_current)
        self.store_time('vdw current' , dt_current)

        return current, current_grp, current_grp_inner


import lib_current
class StressTensor:

    def __init__(self, target_atoms, iatm_to_igrp, bonded_pairs):
        lib_current.bonded.initialize(target_atoms, iatm_to_igrp, bonded_pairs)
        lib_current.nonbonded.initialize(target_atoms, iatm_to_igrp)
        lib_current.kinetic.initialize(target_atoms, iatm_to_igrp)

    def cal_bonded(self, crd, tbfs, volumes, gvolumes=[]):
        """Calculate the current due to bonded potentials."""
        m_bond = lib_current.bonded

        # calculate
        m_bond.cal_bonded(crd, tbfs, volumes, gvolumes)

        # get current by copy
        cur_atm = m_bond.current.copy()
        cur_inn = m_bond.current_inn.copy()
        cur_out = m_bond.current_out.copy()
        return cur_atm, cur_inn, cur_out

    def cal_nonbonded(self, crd, gen_tbfs, volumes, table, gvolumes=[]):
        """Calculate the current due to nonbonded potentials."""
        t0 = time.time()

        m_non = lib_current.nonbonded

        # initialize
        m_non.init_cal(crd, volumes, gvolumes)

        # calculate
        t_total = time.time() - t0
        for t, tbfs in zip(table, gen_tbfs):
            t1 = time.time()
            m_non.cal_nonbonded(tbfs, t)
            t_total += time.time() - t1

        t2 = time.time()

        # get current by copy
        cur_atm = m_non.current.copy()
        cur_inn = m_non.current_inn.copy()
        cur_out = m_non.current_out.copy()

        t_total += time.time() - t2
        self.dt = t_total

        return cur_atm, cur_inn, cur_out

    def cal_kinetic(self, vel, masses, volumes, gvolumes=[]):
        """Calculate the current due to kinetic terms."""
        m_kin = lib_current.kinetic

        # convert the unit from g*[A/fs]^2 to [kcal/mol]
        s = 10.0**4 / 4.184

        # calculate
        m_kin.cal_kinetic(vel, masses, volumes, gvolumes)

        # get current
        return m_kin.current*s, m_kin.current_inn*s, m_kin.current_out*s

    def cal_bonded_old(self, crd, tbfs, volumes):
        """Calculate the current due to bonded potentials."""
        return lib_current.cal_bonded(crd, tbfs, volumes)

    def cal_nonbonded_old(self, crd, gen_tbfs, volumes, table):
        """Calculate the current due to nonbonded potentials."""

        # initialize
        lib_current.nonbond.initialize(crd, volumes)

        for t, tbfs in zip(table, gen_tbfs):
            lib_current.nonbond.cal_nonbonded(tbfs, t)

        return lib_current.nonbond.current.copy()

    def cal_kinetic_old(self, vel, volumes, masses):
        ps = [m*v/vol for v, vol, m in zip(vel, volumes, masses)]
        return utility.tensor(numpy.array(ps), vel)



################################################################################

# TODO: Theory is yet not created ...
class EnergyCurrentCalculator(base.CurrentCalculator):

    def cal_bonded(self):
        """Calculate the energy current for the bonded potential terms."""
        
        # calculate two-body force of bonded interactions
        self.__tbf.cal_bonded()

        # get two-body forces
        bonded_tbfs = self.__tbf.get_bonded_tb_forces()
        bonded14_tbfs = self.__tbf.get_bonded14_tb_forces()

        # initialize
        current_pot = numpy.zeros([self.__natom+1, 3])

        # bonded two-body force
        for iatm, jatm, f_ij in bonded_tbfs.items():
            v_ij = 0.5 * (vel[iatm] + vel[jatm])
            f_ij_v_ij = f_ij * v_ij
            current_pot[iatm] += f_ij_v_ij / volume_ij
            current_pot[jatm] += f_ij_v_ij / volume_ji

        # bonded14 two-body force
        for iatm, jatm, f_ij in bonded14_tbfs.items():
            v_ij = 0.5 * (vel[iatm] + vel[jatm])
            f_ij_v_ij = f_ij * v_ij
            current_pot[iatm] += f_ij_v_ij / volume_ij
            current_pot[jatm] += f_ij_v_ij / volume_ji

        return current_pot

    #TODO: volume calculation
    def _cal_nonbond(self):
        current_pot = numpy.zeros([self.__natom+1, 3, 3])

        for iatm, jatm_ptr, f_ijs in gen_tbforces_each_atom:

            for j, f_ij in enumerate(f_ijs):
                jatm = jatm_ptr + j

                v_ij = 0.5 * (vel[iatm] + vel[jatm])
                f_ij_v_ij = f_ij * v_ij
                current_pot[iatm] += f_ij_v_ij / volume_ij
                current_pot[jatm] += f_ij_v_ij / volume_ji

        return current_pot
                

if __name__ == '__main__':
    cal = StressCurrentCalculator()
    # table = cal.make_interaction_table(natom=5)
