from __future__ import print_function

import os, sys
import time
import numpy

# curp modules
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)
import utility
import clog as logger

import base
################################################################################
class EnergyFluxCalculator(base.FluxCalculator):

    """
    This class calculates the energy flux.
    """

    def __init__(self):
        base.FluxCalculator.__init__(self)

    def prepare(self, *args, **kwds):
        """Prepare the calculations of all snapshots."""
        base.FluxCalculator.prepare(self, *args, **kwds)

        ftype = self.get_setting().curp.flux_grain
        if ftype == 'atom':
            flag_atom  = True
            flag_group = False
        elif ftype == 'group':
            flag_atom  = False
            flag_group = True
        elif ftype == 'both':
            flag_atom  = True
            flag_group = True
        else:
            pass

        self.fcal = EnergyFlux( self.get_target_atoms(),
                self.get_iatm_to_igrp(), self.get_bonded_pairs(),
                flag_atom, flag_group)

    def cal_bonded(self, crd, vel, bond_type):
        """Calculate the energy flux for the bonded term."""

        tbcal = self.get_tbforce()

        # calculate two-body force of bonded interactions
        t0 = time.time()
        tbfs_info = tbcal.cal_bonded(bond_type)
        tbfs      = tbfs_info['tbforces']


        # calculate current
        t1 = time.time()
        flux_atm, flux_grp = self.fcal.cal_bonded(vel, tbfs)

        # log the elasped times
        t2 = time.time()
        self.store_time(bond_type+' pairwise', t1-t0)
        self.store_time(bond_type+' flux' , t2-t1)

        return flux_atm, flux_grp

    def cal_coulomb(self, crd, vel):
        """Calculate the energy flux for the coulomb term."""

        table     = self.get_interact_table()
        type_func = self.get_tbforce().cal_coulomb
        cutoff    = self.get_setting().curp.coulomb_cutoff_length

        t0 = time.time()
        gen_tbfs = ( type_func(t)['tbforces'] for t in table )
        flux_atm, flux_grp = self.fcal.cal_nonbonded( vel, gen_tbfs, table )

        t1 = time.time()
        dt_flux = self.fcal.dt

        self.store_time('coulomb pairwise' , t1 - t0 - dt_flux)
        self.store_time('coulomb flux' , dt_flux)

        return flux_atm, flux_grp

    def cal_vdw(self, crd, vel):
        """Calculate the energy flux for the vdW term."""

        table     = self.get_interact_table()
        type_func = self.get_tbforce().cal_vdw
        cutoff    = self.get_setting().curp.vdw_cutoff_length

        t0 = time.time()
        gen_tbfs = ( type_func(t)['tbforces'] for t in table )
        flux_atm, flux_grp = self.fcal.cal_nonbonded( vel, gen_tbfs, table )

        t1 = time.time()
        dt_flux = self.fcal.dt

        self.store_time('vdw pairwise' , t1 - t0 - dt_flux)
        self.store_time('vdw flux' , dt_flux)

        return flux_atm, flux_grp


################################################################################
import lib_flux
class EnergyFlux:

    def __init__(self, target_atoms, iatm_to_igrp, bonded_pairs,
                       flag_atm=True, flag_grp=True):
        self.__flag_atm = flag_atm
        self.__flag_grp = flag_grp

        lib_flux.bonded.initialize( target_atoms, iatm_to_igrp,
                bonded_pairs, flag_atm, flag_grp)
        lib_flux.nonbonded.initialize( target_atoms, iatm_to_igrp,
                flag_atm, flag_grp)

    def cal_bonded(self, vel, tbfs):
        """Calculate the flux due to bonded potentials."""
        m_bond = lib_flux.bonded

        # calculate
        m_bond.cal_bonded(vel, tbfs)

        # get current by copy
        # flux_atm = m_bond.flux_atm if self.__flag_atm else None
        # flux_grp = m_bond.flux_grp if self.__flag_grp else None
        flux_atm = m_bond.flux_atm.copy() if self.__flag_atm else None
        flux_grp = m_bond.flux_grp.copy() if self.__flag_grp else None
        return flux_atm, flux_grp

    def cal_nonbonded(self, vel, gen_tbfs, table):
        """Calculate the flux due to nonbonded potentials."""
        t0 = time.time()
        m_non = lib_flux.nonbonded

        # initialize
        m_non.init_cal(vel)

        # calculate
        t_total = time.time() - t0
        tt = 0.0
        for t, tbfs in zip(table, gen_tbfs):
            t1 = time.time()
            m_non.cal_nonbonded(tbfs, t)
            t_total += time.time() - t1
            tt += time.time() - t1
            # print('flux body loop: ', time.time()-t1)

        # print('flux body:', tt)
        t2 = time.time()

        # get current by copy
        # flux_atm = m_non.flux_atm if self.__flag_atm else None
        # flux_grp = m_non.flux_grp if self.__flag_grp else None
        flux_atm = m_non.flux_atm.copy() if self.__flag_atm else None
        flux_grp = m_non.flux_grp.copy() if self.__flag_grp else None

        t_total += time.time() - t2
        self.dt = t_total

        return flux_atm, flux_grp

################################################################################
################################## HEAT FLUX ###################################

class HeatFluxCalculator(base.FluxCalculator):

    """
    This class calculates the heat flux.
    """

    def __init__(self):
        base.FluxCalculator.__init__(self)

    def prepare(self, *args, **kwds):
        """Prepare the calculations of all snapshots."""
        base.FluxCalculator.prepare(self, *args, **kwds)

        ftype = self.get_setting().curp.flux_grain
        if ftype == 'atom':
            flag_atom  = True
            flag_group = False
        elif ftype == 'group':
            flag_atom  = False
            flag_group = True
        elif ftype == 'both':
            flag_atom  = True
            flag_group = True
        else:
            pass

        self.fcal = HeatFlux( self.get_target_atoms(),
                self.get_iatm_to_igrp(), self.get_bonded_pairs(),
                flag_atom, flag_group)

    def cal_bonded(self, crd, vel, bond_type):
        """Calculate the energy flux for the bonded term."""

        tbcal = self.get_tbforce()

        # calculate two-body force of bonded interactions
        t0 = time.time()
        tbfs_info = tbcal.cal_bonded(bond_type)
        tbfs      = tbfs_info['tbforces']
        displacement = tbfs_info['displacement']

        # calculate current
        t1 = time.time()
        flux_atm, flux_grp = self.fcal.cal_bonded(vel, tbfs, displacement)

        # log the elasped times
        t2 = time.time()
        self.store_time(bond_type+' pairwise', t1-t0)
        self.store_time(bond_type+' flux' , t2-t1)

        return flux_atm, flux_grp

    def cal_coulomb(self, crd, vel):
        """Calculate the energy flux for the coulomb term."""

        table     = self.get_interact_table()
        type_func = self.get_tbforce().cal_coulomb
        cutoff    = self.get_setting().curp.coulomb_cutoff_length

        t0 = time.time()
        gen_tbfs = (type_func(t)['tbforces'] for t in table)
        gen_displacement = ( type_func(t)['displacement'] for t in table )

        flux_atm, flux_grp = self.fcal.cal_nonbonded( vel, gen_tbfs, table,
                                                      gen_displacement )

        t1 = time.time()
        dt_flux = self.fcal.dt

        self.store_time('coulomb pairwise' , t1 - t0 - dt_flux)
        self.store_time('coulomb flux' , dt_flux)

        return flux_atm, flux_grp

    def cal_vdw(self, crd, vel):
        """Calculate the energy flux for the vdW term."""

        table     = self.get_interact_table()
        type_func = self.get_tbforce().cal_vdw
        cutoff    = self.get_setting().curp.vdw_cutoff_length

        t0 = time.time()
        gen_tbfs = ( type_func(t)['tbforces'] for t in table )
        gen_displacement = ( type_func(t)['displacement'] for t in table )

        flux_atm, flux_grp = self.fcal.cal_nonbonded( vel, gen_tbfs, table,
                                                      gen_displacement )

        t1 = time.time()
        dt_flux = self.fcal.dt

        self.store_time('vdw pairwise' , t1 - t0 - dt_flux)
        self.store_time('vdw flux' , dt_flux)

        return flux_atm, flux_grp


################################################################################
import lib_hflux
class HeatFlux:

    def __init__(self, target_atoms, iatm_to_igrp, bonded_pairs,
                       flag_atm=True, flag_grp=True):
        self.__flag_atm = flag_atm
        self.__flag_grp = flag_grp

        lib_hflux.bonded.initialize( target_atoms, iatm_to_igrp,
                bonded_pairs, flag_atm, flag_grp)
        lib_hflux.nonbonded.initialize( target_atoms, iatm_to_igrp,
                flag_atm, flag_grp)

    def cal_bonded(self, vel, tbfs, displ):
        """Calculate the flux due to bonded potentials."""
        m_bond = lib_hflux.bonded
        # calculate
        m_bond.cal_bonded(vel, tbfs, displ)

        # get current by copy
        # flux_atm = m_bond.flux_atm if self.__flag_atm else None
        # flux_grp = m_bond.flux_grp if self.__flag_grp else None
        hflux_atm = m_bond.hflux_atm.copy() if self.__flag_atm else None
        hflux_grp = m_bond.hflux_grp.copy() if self.__flag_grp else None
        return hflux_atm, hflux_grp

    def cal_nonbonded(self, vel, gen_tbfs, table, gen_displ):
        """Calculate the flux due to nonbonded potentials."""
        t0 = time.time()
        m_non = lib_hflux.nonbonded

        # initialize
        m_non.init_cal(vel)

        # calculate
        t_total = time.time() - t0
        tt = 0.0
        for t, tbfs, displ in zip(table, gen_tbfs, gen_displ):
            t1 = time.time()
            m_non.cal_nonbonded(tbfs, t, displ)
            t_total += time.time() - t1
            tt += time.time() - t1
            # print('flux body loop: ', time.time()-t1)

        # print('flux body:', tt)
        t2 = time.time()

        # get current by copy
        # flux_atm = m_non.flux_atm if self.__flag_atm else None
        # flux_grp = m_non.flux_grp if self.__flag_grp else None
        hflux_atm = m_non.hflux_atm.copy() if self.__flag_atm else None
        hflux_grp = m_non.hflux_grp.copy() if self.__flag_grp else None

        t_total += time.time() - t2
        self.dt = t_total

        return hflux_atm, hflux_grp




################################################################################
#TODO
class StressFluxCalculator(base.FluxCalculator):

    #TODO: volume calculation
    def cal_bonded_flux(self):
        """Calculate the stress flux for the bonded terms."""

        # calculate two-body force of bonded interactions
        self.__tbf.cal_bonded()

        # get two-body forces
        bonded_tbfs = self.__tbf.get_bonded_tb_forces()
        bonded14_tbfs = self.__tbf.get_bonded14_tb_forces()

        # initialize
        group_names = self.get_groupnames()
        flux_pot = {}
        for gname_i in group_names:
            for gname_j in group_names:
                ts = numpy.zeros([3,3])
                flux_pot[gname_i, gname_j] = ts

        # bonded two-body force
        for iatm, jatm, f_ij in bonded_tbfs.items():
            gname_i = self.__iatm_to_groups[iatm]
            gname_j = self.__iatm_to_groups[jatm]

            r_ij = crd[iatm] - crd[jatm]
            # f_ij_r_ij = utility.tensor1(f_ij, r_ij)
            flux_pot[gname_i, gname_j] += f_ij_r_ij / volume_ij
            flux_pot[gname_j, gname_i] += f_ij_r_ij / volume_ji

        # bonded14 two-body force
        for iatm, jatm, f_ij in bonded14_tbfs.items():
            gname_i = self.__iatm_to_groups[iatm]
            gname_j = self.__iatm_to_groups[jatm]

            r_ij = crd[iatm] - crd[jatm]
            # f_ij_r_ij = utility.tensor1(f_ij, r_ij)
            flux_pot[gname_i, gname_j] += f_ij_r_ij / volume_ij
            flux_pot[gname_j, gname_i] += f_ij_r_ij / volume_ji

        return flux_pot

    def _cal_nonbond_flux(self):

        # initialize
        group_names = self.get_groupnames()
        flux_pot = {}
        for gname_i in group_names:
            for gname_j in group_names:
                ts = numpy.zeros([3,3])
                flux_pot[gname_i, gname_j] = ts

        for iatm, jatm_ptr, f_ijs in gen_tbforces_each_atom:

            for j, f_ij in enumerate(f_ijs):
                jatm = jatm_ptr + j

                gname_i = self.__iatm_to_groups[iatm]
                gname_j = self.__iatm_to_groups[jatm]

                crd_ij = crd[iatm] - crd[jatm]
                # f_ij_r_ij = utility.tensor1(f_ij, crd_ij)
                flux_pot[gname_i, gname_j] += f_ij_r_ij / volume_ij
                flux_pot[gname_j, gname_i] += f_ij_r_ij / volume_ji

        return flux_pot



if __name__ == '__main__':
    cal = EnergyFluxCalculator()
