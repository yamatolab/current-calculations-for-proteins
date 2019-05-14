from __future__ import print_function

# standard module
import os, sys
import numpy
import time

# curp module
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)
import clog as logger

################################################################################
class TwoBodyForceBase:

    def __init__(self, topology, setting=None):
        self.__tpl = topology
        self.__setting = setting
        self.__natom = self.__tpl.get_natom()
        self.__forces = None
        self.__ptype_to_energy = {}
        self.__ptype_to_forces = {}
        self.__ptype_to_displacement = {}

    def get_pottypes(self):
        return ['bond','angle','torsion','improper',
                'coulomb14','vdw14','coulomb','vdw']
    
    def set_module(self, module):
        self.__mod = module

    def get_module(self):
        return self.__mod

    def get_natom(self):
        return self.__natom

    def setup(self, interact_table, check=False):
        self.__interact_table = interact_table
        max_tbf = self.get_maxpair(interact_table)
        self._setup_init(max_tbf, check)

        self._setup_bond()
        self._setup_angle()
        self._setup_torsion()
        self._setup_improper()

        self._setup_coulomb14()
        self._setup_vdw14()
        self._setup_coulomb()
        self._setup_vdw()

    def cal_force(self, crd):
        # initialize
        self.initialize(crd)

        # calculate the bonded components.
        for ptype in ('bond','angle','torsion','improper','coulomb14','vdw14'):
            self.cal_bonded(ptype)

        # calculate the nonbonded components.
        for t in self.__interact_table:
            self.cal_coulomb(t)
            self.cal_vdw(t)

        return self.__forces

    def _setup_init(self, max_tbf, check=False):
        self.__mod.setup(natom=self.get_natom(), check=check,
                bonded_pairs = self.__tpl.get_bonded_pairs(), max_tbf=max_tbf)

    def _setup_bond(self):
        """Prepare the parameter for the bond calculation."""
        mod = self.__setup_bondtype('bond')
        mod.ibnd_to_itbf = self.__tpl.get_ibnd_to_ipair()

    def _setup_angle(self):
        """Prepare the parameter for the angle calculation."""
        mod = self.__setup_bondtype('angle')
        mod.iang_to_itbf = self.__tpl.get_iang_to_ipair()

    def _setup_torsion(self):
        """Prepare the parameter for the torsion calculation."""
        mod = self.__setup_bondtype('torsion')
        mod.itor_to_itbf = self.__tpl.get_itor_to_ipair()

    def _setup_improper(self):
        """Prepare the parameter for the improper torsion calculation."""
        mod = self.__setup_bondtype('improper')
        mod.itor_to_itbf = self.__tpl.get_iimp_to_ipair()

    def __setup_bondtype(self, btype_name):
        """Prepare the parameter for the calculations without coulomb and vdw.
        """
        info = getattr(self.__tpl, 'get_'+btype_name+'_info')()
        mod_type = getattr(self.__mod, btype_name)
        for key in info.keys():
            value = info[key]
            setattr(mod_type, key, value)

        if btype_name == 'torsion':
            logger.info('The number of torsions :', len(mod_type.four_atoms))

        if btype_name == 'improper':
            logger.info('The number of impropers :', len(mod_type.four_atoms))

        return mod_type

    def _setup_coulomb(self):
        """Prepare the parameter for the coulomb calculation."""
        coulomb = self.__mod.coulomb
        info = self.__tpl.get_coulomb_info()
        coulomb.charges = info['charges']
        coulomb.cutoff_length = self.__setting.curp.coulomb_cutoff_length

    def _setup_vdw(self):
        """Prepare the parameter for the vdw calculation."""
        vdw = self.__mod.vdw
        info = self.__tpl.get_vdw_info()
        vdw.atom_types = info['atom_types']
        vdw.c6s        = info['c6s']
        vdw.c12s       = info['c12s']
        vdw.cutoff_length = self.__setting.curp.vdw_cutoff_length

    def _setup_coulomb14(self):
        """Prepare the parameter for the coulomb calculation."""
        coulomb14 = self.__mod.coulomb14
        info = self.__tpl.get_coulomb_info()
        coulomb14.charges = info['charges']
        # return self.__tpl.get_i14_to_ipair()
        i14_to_itbf = self.__tpl.get_i14_to_ipair()

        # if len(self.__tpl.get_i14_to_ipair())==1:
            # i14_to_itbf = numpy.array([])
        # else:
            # i14_to_itbf = self.__tpl.get_i14_to_ipair()

        # coulomb14.setup(info['charges'], i14_to_itbf)
        coulomb14.i14_to_itbf = self.__tpl.get_i14_to_ipair()

    def _setup_vdw14(self):
        """Prepare the parameter for the vdw calculation."""
        vdw14 = self.__mod.vdw14
        info = self.__tpl.get_vdw_info()
        vdw14.atom_types = info['atom_types']
        vdw14.c6s        = info['c6s']
        vdw14.c12s       = info['c12s']

        # if len(self.__tpl.get_i14_to_ipair())==1:
            # i14_to_itbf = numpy.array([])
        # else:
        i14_to_itbf = self.__tpl.get_i14_to_ipair()

        # vdw14.setup(info['atom_types'],info['c6s'],info['c12s'], i14_to_itbf)
        # print(self.__tpl.get_i14_to_ipair())
        vdw14.i14_to_itbf = self.__tpl.get_i14_to_ipair()
        # print(vdw14.i14_to_itbf )

        logger.info( 'The number of 1-4 interactions', len(i14_to_itbf))

    def initialize(self, crd):
        self.__mod.initialize(crd)
        self.__forces   = numpy.zeros( [self.__natom, 3] )
        self.__ptype_to_energy = {}
        self.__ptype_to_forces = {}

    def cal_bonded(self, bond_type):
        """Calculate the pairwise forces using the bonded type modules.
        Types are bond, angle, torsion, improper, coulomb14 and vdw14.
        """
        mod = getattr(self.__mod, bond_type)
        mod.calculate()
        self.__forces += mod.forces

        # store energy and forces
        self.__ptype_to_energy[bond_type] = mod.energy
        self.__ptype_to_forces[bond_type] = mod.forces
        self.__ptype_to_displacement[bond_type] = mod.displacement

        #DEBUG
        # print('** {} forces **'.format(bond_type))
        # for iatm_1, f in enumerate(mod.forces):
            # print("{:5>} {:15.7f}{:15.7f}{:15.7f}".format(
                # iatm_1+1, f[0],f[1],f[2]))
        # print()

        # if bond_type == 'bond':
            # print('** {} pairwise forces **'.format(bond_type))
            # for ibnd_1, itbf in enumerate(mod.ibnd_to_itbf):
                # tbf = mod.tbforces[itbf-1]
                # iatm, jatm = mod.two_atoms[ibnd_1]
                # print("{:5>} {:5>} {:5>} {:15.7f}{:15.7f}{:15.7f}".format(
                    # ibnd_1+1, iatm,jatm, tbf[0],tbf[1],tbf[2]))
            # print()

        #DEBUG

        return dict(energy = mod.energy,
                    forces = mod.forces,
                    tbforces = mod.tbforces,
                    displacement = mod.displacement)

    def cal_coulomb(self, table):
        return self._cal_nonbond(table, 'coulomb')

    def cal_vdw(self, table):
        return self._cal_nonbond(table, 'vdw')

    def _cal_nonbond(self, table, pottype):
        """Calculate the pairwise forces using the bonded type modules.
        Types are coulomb and vdw.
        """

        mod = getattr(self.__mod, pottype)
        mod.calculate(table)
        energy = mod.energy.copy()
        forces = mod.forces.copy()
        displacement = mod.displacement

        # store energy, forces and distance
        if pottype not in self.__ptype_to_energy:
            self.__ptype_to_energy[pottype] = energy
        else:
            self.__ptype_to_energy[pottype] += energy

        self.__forces += forces
        if pottype not in self.__ptype_to_forces:
            self.__ptype_to_forces[pottype] = forces
        else:
            self.__ptype_to_forces[pottype] += forces

        self.__ptype_to_displacement[pottype] = displacement

        # get pairwise forces
        tbforces = mod.tbforces.copy()

        return dict(energy = mod.energy,
                    forces = mod.forces,
                    tbforces = tbforces,
                    displacement = displacement)

    def get_forces(self, pottype):
        """Return the calculated force with givin potential type."""
        return self.__ptype_to_forces[pottype]

    def get_energy(self, pottype):
        """Return the calculated energy."""
        return self.__ptype_to_energy[pottype]

    def output_energy(self):
        """Output the energy."""
        logger.info_cycle('    ** Output energy **')
        for pottype in self.get_pottypes():
            energy = self.get_energy(pottype)
            logger.info_cycle('    {:>10} : {:>}'.format(pottype, energy))
        logger.info_cycle()

    def output_force(self):
        logger.debug_cycle('    ** Output force **')
        for pottype in self.get_pottypes():
            logger.debug_cycle('    [ {} force ] '.format(pottype))

            force = self.get_forces(pottype)
            for iatm_1, f in enumerate(force):
                logger.debug_cycle(
                    '    {:>8} : {:>12.8f} {:>12.8f} {:>12.8f}'.format(
                    iatm_1+1, f[0], f[1], f[2]))

        logger.debug_cycle('    [ Total force ] '.format(pottype))
        for iatm_1, f in enumerate(self.__forces):
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
        caltype_mod = getattr(self.get_module(), pot_type)
        caltype_mod.print_tbforce()
        print()

    # def write_bond(self):
    #     fmt = "{iatm} {jatm} {x} {y} {z}"
    #     tbfs = self.__mod.bond.tbforces
    #     two_atoms = self.__tpl.get_bond_info()['two_atoms']
    #     for two, tbf in zip(two_atoms, tbfs):
    #         iatm, jatm = two
    #         print(fmt.format(
    #             iatm=iatm, jatm=jatm, x=tbf[0], y=tbf[1],z=tbf[2]))

    # def write_energy(self):
    #     pass

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
                print(iatm, jatm, tbfs[itbf-1])

        print()

    def get_maxpair(self, interact_table):
        maxpair = 0
        for t in interact_table:
            npair = 0
            for iatm, jatm_beg, jatm_end in t:
                npair += jatm_end - jatm_beg + 1

            maxpair = max(maxpair, npair)

        return maxpair

class TwoBodyForce(TwoBodyForceBase):

    def __init__(self, topology, setting=None):
        TwoBodyForceBase.__init__(self, topology, setting)
        import lib_amberbase
        self.set_module(lib_amberbase)


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

