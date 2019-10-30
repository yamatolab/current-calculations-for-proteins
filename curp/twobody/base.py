from __future__ import print_function, division, absolute_import

# standard module
import os
import sys
import numpy

# curp module
from .. import clog as logger
from . import lib_base

########################################################################
class TwoBodyForce:
    module = lib_base
    def __init__(self, topology, setting=None):
        self.tpl = topology
        self._setting = setting
        self._forces = None
        self._ptype_to_energy = {}
        self._ptype_to_forces = {}
        self._ptype_to_displacement = {}

        self._natom = self.tpl.natom
        self._pottypes = self.tpl.decomp_list

    @property
    def pottypes(self):
        # Potential energy types
        return self._pottypes

    @property
    def natom(self):
        return self._natom

    def setup(self, interact_table, check=False):
        """Initialize fortran module"""
        self._interact_table = interact_table
        max_tbf = self.get_maxpair(interact_table)
        self.module.setup(natom=self.natom, check=check,
                          bonded_pairs=self.tpl.bonded_pairs, max_tbf=max_tbf)

        self._setup_bonded()
        logger.info('The number of impropers :', len(self.tpl.impropers))

        self._setup_coulomb()
        self._setup_vdw()

    def get_maxpair(self, interact_table):
        maxpair = 0
        for table in interact_table:
            npair = 0
            for iatm, jatm_beg, jatm_end in table:
                npair += jatm_end - jatm_beg + 1

            maxpair = max(maxpair, npair)

        return maxpair

    def cal_force(self, crd):
        """Calculate forces"""
        # Initialize
        self.module.initialize(crd)
        self._forces = numpy.zeros([self._natom, 3])
        self._ptype_to_energy = {}
        self._ptype_to_forces = {}

        # Calculate the bonded components.
        for ptype in self.tpl.decomp_list('bonded+'):
            self.cal_bonded(ptype)

        # Calculate the nonbonded components.
        for table in self._interact_table:
            self.cal_nonbond(table, 'coulomb')
            self.cal_nonbond(table, 'vdw')

        return self._forces

    def _setup_bonded(self):
        """Prepare the parameter for the calculations without coulomb and vdw."""
        bonded_inter = self.tpl.bonded_inter
        for btype, interactions in bonded_inter.items():
            int_mod = getattr(self.module, btype)
            # Get the name of the constants and lists of values.
            for cst_name, values in interactions.ff_cst.items():
                setattr(int_mod, cst_name, values)
            int_mod.to_ipair = interactions.to_ipair
            int_mod.ids = interactions.ids

    def _setup_coulomb(self):
        """Prepare the parameter for the coulomb calculation."""
        coulomb = self.module.coulomb
        coulomb.charges = self.tpl.atoms['charges']
        coulomb.cutoff_length = self._setting.curp.coulomb_cutoff_length

    def _setup_vdw(self):
        """Prepare the parameter for the vdw calculation."""
        vdw = self.module.vdw
        vdw.vdw_radii = self.tpl.atoms['vdw_radii']
        vdw.epsilons = self.tpl.atoms['epsilons']
        vdw.cutoff_length = self._setting.curp.vdw_cutoff_length

    def cal_bonded(self, bond_type):
        """Calculate the pairwise forces using the bonded type modules.
        Types are bond, angle, dihedral and improper.
        """
        mod = getattr(self.module, bond_type)
        mod.calculate()
        self._forces += mod.forces

        # Store energy and forces
        self._ptype_to_energy[bond_type] = mod.energy
        self._ptype_to_forces[bond_type] = mod.forces
        self._ptype_to_displacement[bond_type] = mod.displacement

        return dict(energy = mod.energy,
                    forces = mod.forces,
                    tbforces = mod.tbforces,
                    displacement = mod.displacement)

    def _cal_nonbond(self, table, pottype):
        """Calculate the pairwise forces using the fortran module.

        pottype either coulomb or vdw.
        """

        mod = getattr(self.module, pottype)
        mod.calculate(table)
        energy = mod.energy.copy()
        forces = mod.forces.copy()
        displacement = mod.displacement

        # store energy, forces and distance
        if pottype not in self._ptype_to_energy:
            self._ptype_to_energy[pottype] = energy
        else:
            self._ptype_to_energy[pottype] += energy

        self._forces += forces
        if pottype not in self._ptype_to_forces:
            self._ptype_to_forces[pottype] = forces
        else:
            self._ptype_to_forces[pottype] += forces

        self._ptype_to_displacement[pottype] = displacement

        # get pairwise forces
        tbforces = mod.tbforces.copy()

        return dict(energy = mod.energy,
                    forces = mod.forces,
                    tbforces = tbforces,
                    displacement = displacement)

    @property
    def forces(self):
        """Return the calculated force for each potential type."""
        return self._ptype_to_forces

    @property
    def energies(self):
        """Return the calculated energy."""
        return self._ptype_to_energy

    def output_energy(self):
        """Output the energy."""
        logger.info_cycle('    ** Output energy **')
        for pottype in self.pottypes:
            energy = self.energies[pottype]
            logger.info_cycle('    {:>10} : {:>}'.format(pottype, energy))
        logger.info_cycle()

    def output_force(self):
        logger.debug_cycle('    ** Output force **')
        for pottype in self.pottypes:
            logger.debug_cycle('    [ {} force ] '.format(pottype))

            force = self.forces[pottype]
            for iatm_1, f in enumerate(force):
                logger.debug_cycle(
                    '    {:>8} : {:>12.8f} {:>12.8f} {:>12.8f}'.format(
                    iatm_1+1, f[0], f[1], f[2]))

        logger.debug_cycle('    [ Total force ] '.format(pottype))
        for iatm_1, f in enumerate(self._forces):
            logger.debug_cycle(
                '    {:>8} : {:>12.8f} {:>12.8f} {:>12.8f}'.format(
                iatm_1+1, f[0], f[1], f[2]))

        logger.debug_cycle()

    def output_bonded(self, results, pot_type):
        print("[{}]".format(pot_type))
        print("== Energy = ")
        print(results['energy'])
        print()

        print("== Forces ==")
        print(results['forces'])
        print()

        print('== Two-body forces ==')
        caltype_mod = getattr(module, pot_type)
        caltype_mod.print_tbforce()
        print()

    def output_nonbonded(self, results, pot_type, table):
        print("[{}]".format(pot_type))
        print("== Energy == ")
        print(results['energy'])
        print()

        print("== Forces ==")
        print(results['forces'])
        print()

        print('== Two-body forces ==')
        tbfs = results['tbforces']
        itbf = 0
        for iatm, jatm_beg, jatm_end in table:
            for jatm in range(jatm_beg, jatm_end+1):
                itbf += 1
                print(iatm, jatm, tbfs[itbf-1], '\n')


if __name__ == '__main__':
    import numpy

    class Setting:
        class Curp:
            coulomb_cutoff_length = 10.0
            vdw_cutoff_length = 10.0

        curp = Curp()


    import os, sys
    topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if topdir not in sys.path:
        sys.path.insert(0, topdir)

    from forcefield.amberbase import ConverterBase

    class DummyTopology(ConverterBase):

        def __init__(self):
            ConverterBase.__init__(self, None)

        def get_mol_info(self): pass
        def convert(self): pass
        def get_residue_info(self): pass
        def get_pbc_info(self): pass

        def get_natom(self):
            return 5

        def get_bond_info(self):
            return dict(
                two_atoms    = [[1,2], [1,3], [4,5]],
                force_consts = [1.0, 1.5, 1.9],
                length_eqs   = [1.0, 1.2, 0.9] )

        def get_angle_info(self):
            return dict(
                three_atoms  = [[1,2,3], [2,3,4], [1,2,5]],
                force_consts = [1.5, 2.0, 3.0],
                theta_eqs    = [90.0, 105.0, 120.0] )

        def get_torsion_info(self):
            return dict(
                four_atoms     = [[1,2,3,4], [2,1,4,5]],
                num_torsions   = [1, 1],
                num_freqs      = [2, 5],
                force_consts   = [1.5, 1.0],
                initial_phases = [180.0, 0.0] )

        def get_improper_info(self):
            return dict(
                four_atoms     = [[1,2,3,4], [2,1,4,5]],
                num_torsions   = [1, 1],
                num_freqs      = [2, 5],
                force_consts   = [1.5, 1.0],
                initial_phases = [180.0, 0.0] )

        def get_coulomb_info(self):
            return dict(
                charges = [0.5, 0.1, 0.3, 0.01, -0.6] )

        def get_vdw_info(self):
            return dict(
                atom_types = [1, 2, 1, 3, 1],
                c6s = [[100.0, 50.0, 35.5],
                       [50.0, 200.0, 10.0],
                       [35.0, 10.0, 300.0]],
                c12s = [[300.0, 91.0, 25.0],
                       [91.0, 200.0, 15.0],
                       [25.0, 15.0, 100.0]] )

        def get_bonded14_pairs(self):
            return [[1,4],[2,5]]

    # setup
    tbcal = TwoBodyForce(DummyTopology(), Setting())
    tbcal.setup(max_tbf=10)

    # # initialize with coordinate
    crd = numpy.array(
        [[1.0, 2.0, 3.0],
         [2.0, 1.0, 1.5],
         [3.0, 3.0, 3.0],
         [0.0, 0.0, 0.0],
         [5.5, 5.1, 4.5]])
    tbcal.initialize(crd)

    # # calculate two-body forces and output its data
    l = ['bond', 'angle', 'torsion', 'improper', 'coulomb14', 'vdw14']
    # l = ['coulomb14', 'vdw14']
    for caltype in l:
        res = tbcal.cal_bonded(caltype)
        tbcal.output_bonded(res, caltype)

    table = [[1, 2, 5], [2, 3, 5], [3,4,5], [4, 5, 5] ]
    for caltype in ['coulomb', 'vdw']:
        res = getattr(tbcal, 'cal_'+caltype)(table)
        tbcal.output_nonbonded(res, caltype, table)

    # crd = numpy.array(
    #     [[2.0, 3.0, 4.1],
    #      [3.0, 2.0, 2.5],
    #      [4.0, 4.0, 4.0],
    #      [1.0, 1.0, 1.0],
    #      [9.0, 9.1, 9.2]])
    # tbcal.initialize(crd, check=True)
    # tbcal.setup()
    # results = tbcal.cal_bonded()
    # output(results)

    # from benchmarker import Benchmarker
    # with Benchmarker(width=20) as bm:

    #     with bm('setup'):
    #         crd = numpy.array(
    #             [[2.0, 3.0, 4.1],
    #              [3.0, 2.0, 2.5],
    #              [4.0, 4.0, 4.0],
    #              [1.0, 1.0, 1.0],
    #              [9.0, 9.1, 9.2]])
    #         tbcal.initialize(crd, check=False)
    #         tbcal.setup()

    #     with bm('bonded'):
    #         results = tbcal.cal_bonded()
    #         output(results)

    #     table = [[1, 2, 5], [2, 3, 3], [2, 5, 5]]
    #     with bm('coulomb'):
    #         results = tbcal.cal_coulomb(table)
    #         output(results)

    #     with bm('vdw'):
    #         results = tbcal.cal_vdw(table)
    #         output(results)

