from __future__ import print_function

import os, sys
import time
import numpy

# curp modules
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path: sys.path.insert(0, topdir)
from  utility import TimeStore
import clog as logger

################################################################################
class CalculatorBase(TimeStore):

    def __init__(self):
        TimeStore.__init__(self)
        self.__base_table = None

    def cal_twobody_force(self):
        pass

    def prepare(self, topology, setting, target_atoms,
            gname_iatoms_pairs, interact_table):

        self.__tpl = topology
        self.__setting = setting

        # decide two-body force calculator.
        import twobody
        TwoBodyCalculator = twobody.get_calculator(setting.curp.potential)
        self.__tbf = TwoBodyCalculator(topology, setting)
        self.__tbf.setup(interact_table, check=False)

        # get the number of atoms
        natom = self.get_tbforce().get_natom()

        # for target atoms
        self.__target_atoms = target_atoms

        # for group
        self.__gname_iatoms_pairs = gname_iatoms_pairs
        self.__gnames = [gname for gname, iatoms in gname_iatoms_pairs]

        # interaction table
        self.__interact_table = [ numpy.array(t) for t in interact_table]

        # make a table to convert the atom index into the group index.
        iatm_to_igrp = numpy.zeros([natom])
        for igrp_1, (gname, iatoms) in enumerate(self.__gname_iatoms_pairs):
            for iatm in iatoms:
                iatm_to_igrp[iatm-1] = igrp_1+1

        self.__iatm_to_igrp = iatm_to_igrp

    def get_tbforce(self):
        return self.__tbf

    def get_topology(self):
        return self.__tpl

    def get_table(self):
        return self.__base_table

    def get_setting(self):
        return self.__setting

    def get_gname_iatoms_pairs(self):
        return self.__gname_iatoms_pairs

    def get_groupnames(self):
        return self.__gnames

    def get_iatm_to_igrp(self):
        return self.__iatm_to_igrp

    def get_target_atoms(self):
        return self.__target_atoms

    def get_bonded_pairs(self):
        return self.__tpl.get_bonded_pairs()

    def get_interact_table(self):
        return self.__interact_table

class CurrentCalculator(CalculatorBase):

    def __init__(self):
        CalculatorBase.__init__(self)
        self.__volumes = None
        self.__gvolumes = None

    def prepare(self, *args, **kwds):
        CalculatorBase.prepare(self, *args, **kwds)

        import volume
        self.__volume_obj = volume.get_volume_calculator( self.get_topology(),
                self.get_setting(), self.get_target_atoms(),
                self.get_gname_iatoms_pairs() )

    def get_volume_obj(self):
        return self.__volume_obj

    def run(self, data):
        """Run the current calculations for all of the components."""

        cstep, (crd, vel, pbc) = data
        logger.debug_cycle('    calculating current values at step {} ...'
                .format(cstep))

        # get twobody force object
        tbcal = self.get_tbforce()
        tbcal.initialize(crd)

        # calculate volume
        t0 = time.time()
        volume  = self.cal_volume(crd)
        gvolume = self.cal_gvolume(len(self.get_groupnames()))
        self.store_time('volume', time.time()-t0)

        # print('[Volume]')
        # print(volume)

        # gather the atomic _current for each potential types.
        key_to_acurs = {} # current for atoms
        key_to_icurs = {} # current for inside region of the groups
        key_to_ocurs = {} # current for outside region of the groups

        # kinetic
        cur_atm, cur_inn, cur_out = self.cal_kinetic(vel)
        key_to_acurs['kinetic'] = cur_atm
        key_to_icurs['kinetic'] = cur_inn
        key_to_ocurs['kinetic'] = cur_out

        # bonded
        btypes = self.get_topology().get_decomp_list('bonded+')
        for btype in btypes:
            # bond type
            cur_atm, cur_inn, cur_out = self.cal_bonded(crd, btype)
            key_to_acurs[btype] = cur_atm
            key_to_icurs[btype] = cur_inn
            key_to_ocurs[btype] = cur_out

        # non-bonded
        cur_atm, cur_inn, cur_out = self.cal_coulomb(crd)
        key_to_acurs['coulomb'] = cur_atm
        key_to_icurs['coulomb'] = cur_inn
        key_to_ocurs['coulomb'] = cur_out

        cur_atm, cur_inn, cur_out = self.cal_vdw(crd)
        key_to_acurs['vdw'] = cur_atm
        key_to_icurs['vdw'] = cur_inn
        key_to_ocurs['vdw'] = cur_out

        # group current = inside current + outside current
        key_to_gcurs = {}
        for key in key_to_acurs.keys():
            key_to_gcurs[key] = key_to_icurs[key] + key_to_ocurs[key]

        return cstep, (key_to_acurs, key_to_gcurs, key_to_icurs, key_to_ocurs)

    def cal_volume(self, crd):
        self.__volumes = self.get_volume_obj().get_volume(crd) 
        return self.__volumes

    def get_volume(self):
        return self.__volumes

    def cal_gvolume(self, ngroup):
        self.__gvolumes = self.get_volume_obj().get_gvolume(ngroup) 
        return self.__gvolumes

    def get_gvolume(self):
        return self.__gvolumes


################################################################################
class FluxCalculator(CalculatorBase):
    def __init__(self):
        CalculatorBase.__init__(self)

    def run(self, data):
        """Run the flux calculations for all of the components."""

        cstep, (crd, vel, pbc) = data
        logger.debug_cycle('    calculating flux values at step {} ...'
                .format(cstep))

        # get twobody force object
        tbcal = self.get_tbforce()
        tbcal.initialize(crd)

        # gather the flux for each potential types.
        key_to_aflux = {} # flux for atoms
        key_to_gflux = {} # flux for group

        btypes = self.get_topology().get_decomp_list('bonded+')
        for btype in btypes:
            # check for the amount of improper torsion and torsion.
            # if btype in ['improper','torsion']: continue
            # bond type
            flux_atm, flux_grp = self.cal_bonded(crd, vel, btype)
            key_to_aflux[btype] = flux_atm
            key_to_gflux[btype] = flux_grp

        # non-bonded
        flux_atm, flux_grp = self.cal_coulomb(crd, vel)
        key_to_aflux['coulomb'] = flux_atm
        key_to_gflux['coulomb'] = flux_grp
        flux_atm, flux_grp = self.cal_vdw(crd, vel)
        key_to_aflux['vdw'] = flux_atm
        key_to_gflux['vdw'] = flux_grp
        # total for atom
        if flux_atm is not None:
            total_atm = numpy.zeros( key_to_aflux['vdw'].shape )

            for flux in key_to_aflux.values():
                total_atm += flux

        else:
            total_atm = None

        # total for group
        if flux_grp is not None:
            total_grp = numpy.zeros( key_to_gflux['vdw'].shape )

            for flux in key_to_gflux.values():
                total_grp += flux

        else:
            total_grp = None

        key_to_aflux['total'] = total_atm
        key_to_gflux['total'] = total_grp

        # write energy
        if self.get_setting().output.output_energy:
            self.get_tbforce().output_energy()
        # write forces
        if logger.is_debug():
            self.get_tbforce().output_force()

        return cstep, (key_to_aflux, key_to_gflux)

