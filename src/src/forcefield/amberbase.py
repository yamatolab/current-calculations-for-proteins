from __future__ import print_function

import os, sys
from abc import abstractmethod, abstractproperty, ABCMeta
import numpy

# curp package
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)
import table.interact_table as it
import clog as logger


class ConverterPrintable(object):

    def print_atom(self):
        logger.info('*** atom ***')
        info = self.get_atom_info()
        atom = zip(info['names'], info['elems'], 
                info['masses'], info['vdw_radii'])
        for iatm_1, (name, elem, mass, radius) in enumerate(atom):
            logger.info('{:>7d}  {:<4s} {:<4s} {:12.7f} {:12.7f}'.format(
                iatm_1+1, name, elem, mass, radius))
        logger.info()

    def print_residue(self):
        logger.info('*** residue ***')
        info = self.get_residue_info()
        res = zip(info['ids'], info['names'])
        for iatm_1, (rid, rname) in enumerate(res):
            logger.info('{:>7d} {:>5d} {:5s}'.format(iatm_1+1, rid, rname))
        logger.info()

    def print_bond(self):
        logger.info('*** bond ***')
        info = self.get_bond_info()
        bond = zip(info['two_atoms'], info['force_consts'], info['length_eqs'])
        for ibnd_1, (iatoms, force, eq) in enumerate(bond):
            logger.info('{:>5d} {:>15}  {:12.7f}  {:12.7f}'.format(
                ibnd_1+1, iatoms, force, eq))
        logger.info()

    def print_angle(self):
        logger.info('*** angle ***')
        info = self.get_angle_info()
        angle = zip(info['three_atoms'], info['force_consts'],info['theta_eqs'])
        for iang_1, (iatoms, force, eq) in enumerate(angle):
            logger.info('{:>5d} {:>20}  {:12.7f}  {:12.7f}'.format(
                iang_1+1, iatoms, force, eq))
        logger.info()

    def print_torsion(self):
        logger.info('*** torsion ***')
        info = self.get_torsion_info()
        self._print_torsion(info)

    def print_improper(self):
        logger.info('*** improper ***')
        info = self.get_improper_info()
        self._print_torsion(info)

    def _print_torsion(self, info):
        torsion = zip(info['four_atoms'], info['num_torsions'],
                info['num_freqs'], info['force_consts'], info['initial_phases'])
        for itor_1, (iatoms, ntor, freq, force, phase) in enumerate(torsion):
            logger.info('{:>5d} {:>30}  {} {} {:12.7f}  {:12.7f}'.format(
                itor_1+1, iatoms, ntor, freq, force, phase ))
        logger.info()

    def print_coulomb(self):
        logger.info('*** coulomb ***')
        info = self.get_coulomb_info()
        for iatm_1, charge in enumerate(info['charges']):
            logger.info('{:>7d}  {:12.7f}'.format(iatm_1+1, charge))
        logger.info()

    def print_vdw(self):
        logger.info('*** vdw ***')
        logger.info('{:>5} {:>5}{:>17}{:>17}'
                .format('itype','jtype','c6','c12'))
        info = self.get_vdw_info()
        ntype, tmp = info['c6s'].shape
        c6s = info['c6s']
        c12s = info['c12s']
        for itype_1 in range(ntype):
            for jtype_1 in range(ntype):
                logger.info('{:>5} {:>5}  {:15.9E}  {:15.9E}'.format(itype_1+1,
                    jtype_1+1, c6s[itype_1, jtype_1], c12s[itype_1, jtype_1]))
        logger.info('[atom_types]')
        for iatm_1, itype in enumerate(info['atom_types']):
            logger.info('{:>7}  {:>3}'.format(iatm_1+1, itype))
        logger.info()

    def print_bonded_pairs(self):
        bonded_pairs = self.get_bonded_pairs()
        logger.info('*** bonded pairs ***')
        self._print_pairs(self.get_bonded_pairs())

    def print_bonded14_pairs(self):
        logger.info('*** 14 bonded pairs ***')
        self._print_pairs(self.get_bonded14_pairs())

    def print_nonbonded_table(self):
        logger.info('*** non-bonded interaction table ***')

        interact_table = self.get_nonbonded_table()

        logger.info( '{:>10}{:>10}{:>10}'.format('iatm', 'jatm_beg', 'jatm_end') )
        for iatm, jatm_beg, jatm_end in interact_table:
            logger.info( '{:>10}{:>10}{:>10}'.format(iatm,jatm_beg,jatm_end) )

        logger.info()

    def print_ibnd_to_ipair(self): self._print_idx_to_ipair('bnd')
    def print_iang_to_ipair(self): self._print_idx_to_ipair('ang')
    def print_itor_to_ipair(self): self._print_idx_to_ipair('tor')
    def print_iimp_to_ipair(self): self._print_idx_to_ipair('imp')
    def print_i14_to_ipair(self):  self._print_idx_to_ipair('14')

    def _print_pairs(self, pairs):
        npair = len(pairs)
        logger.info(' The number of pairs = {}'.format(npair))
        logger.info('{:>7} {:>7} {:>7}'.format('ipair', 'iatm', 'jatm'))
        for ipair_1, (iatm, jatm) in enumerate(pairs):
            logger.info('{:>7d} {:>7d} {:>7d}'.format(ipair_1+1, iatm, jatm))
        logger.info()

    def _print_idx_to_ipair(self, bondtype):
        
        logger.info('*** i{} => ipair ***'.format(bondtype))
        get = getattr(self, 'get_i'+bondtype+'_to_ipair')
        idx_to_ipair = get()

        n = len(idx_to_ipair)
        logger.info(' The number of indices = {}'.format(n))

        if len(numpy.shape(idx_to_ipair)) == 1:
            ncol = 1
        else:
            ncol = len(idx_to_ipair[0])
        fmt = '{:>7}' + '{:>' + str(ncol*8) + '}'

        logger.info(fmt.format('index', 'ipairs'))
        for idx_1, ipairs in enumerate(idx_to_ipair):
            logger.info(fmt.format(idx_1+1, ipairs))
        logger.info()


import lib_bonded_pair
class TableMaker:

    def __init__(self):

        self._bonded_pairs    = None
        self._bonded14_pairs  = None
        self._nonbonded_table = None

        self._ibnd_to_ipair = None
        self._iang_to_ipair = None
        self._itor_to_ipair = None
        self._iimp_to_ipair = None
        self._i14_to_ipair  = None

        self._target_atoms = None

    def get_decomp_list(self, dtype='all'):
        """Get the list that decompose all potential.
        dtype = 'all'(defaulst), 'bonded', 'bonded14', 'bonded+',
        'nonbonded or 'nonbonded+(nonbonded+14)'
        """
        bonded_list    = ['bond', 'angle', 'torsion', 'improper']
        bonded14_list  = ['coulomb14', 'vdw14']
        nonbonded_list = ['coulomb', 'vdw']

        if   dtype == 'bonded':     return bonded_list
        elif dtype == 'bonded14':   return bonded14_list
        elif dtype == 'bonded+':    return bonded_list + bonded14_list
        elif dtype == 'nonbonded':  return nonbonded_list
        elif dtype == 'nonbonded+': return nonbonded_list + bonded14_list
        else: return bonded_list + bonded14_list + nonbonded_list

    def get_bonded_pairs(self):
        if self._bonded_pairs is None:
            self._bonded_pairs = self._make_bonded_pairs()
        return self._bonded_pairs

    def get_bonded14_pairs(self):
        if self._bonded14_pairs is None:
            self._bonded14_pairs = self._make_bonded14_pairs()
        return self._bonded14_pairs

    def get_nonbonded_table(self):
        if self._nonbonded_table is None:
            self._nonbonded_table = self._make_nonbonded_table()
        return self._nonbonded_table

    def get_ibnd_to_ipair(self): return self._get_idx_to_ipair('bnd')
    def get_iang_to_ipair(self): return self._get_idx_to_ipair('ang')
    def get_itor_to_ipair(self): return self._get_idx_to_ipair('tor')
    def get_iimp_to_ipair(self): return self._get_idx_to_ipair('imp')

    def get_i14_to_ipair(self):
        return self._get_idx_to_ipair('14')
        # if len(self._get_idx_to_ipair('14'))==1:
            # return numpy.array([])
        # else:
            # return self._get_idx_to_ipair('14')

    #TODO
    def apply_target_atoms(self, target_atoms):
        self.__target_atoms = target_atoms

        # bond
        info = self.get_bond_info()
        iatoms_list = info['two_atoms']
        self._apply_target_to_bondtype(iatoms_list, info)

        # angle
        info = self.get_angle_info()
        iatoms_list = info['three_atoms']
        self._apply_target_to_bondtype(iatoms_list, info)

        # torsion
        info = self.get_torsion_info()
        iatoms_list = info['four_atoms']
        self._apply_target_to_bondtype(iatoms_list, info)

        # improper
        info = self.get_improper_info()
        iatoms_list = info['four_atoms']
        self._apply_target_to_bondtype(iatoms_list, info)

        # bonded14
        self._make_bonded14_pairs()

        # bonded
        self._make_bonded_pairs()

    #TODO
    def _make_nonbonded_table(self):
        import time

        bonded_pairs = numpy.array(self.get_bonded_pairs())

        natom = self.get_natom()
        int_table = it.InteractionTable(self.get_natom())

        import lib_nonbond
        mod_tab = lib_nonbond.without_bonded
        mod_tab.setup(bonded_pairs, natom)

        int_table = numpy.array( list(int_table) )
        nonbonded_table, ntable = mod_tab.get_nonbonded_table( int_table )
        new_table = it.InteractionTable(
                base_table=nonbonded_table[:ntable].tolist())

        return new_table

    def _make_bonded_pairs(self):
        two_atoms   = self.get_bond_info()['two_atoms']
        three_atoms = self.get_angle_info()['three_atoms']
        four_atoms  = self.get_torsion_info()['four_atoms']
        four_atoms_imp  = self.get_improper_info()['four_atoms']

        import itertools as it
        pairs = [] # ipair => (iatm, jatm)
        for two in two_atoms:
            iatm, jatm = min(two), max(two)
            pairs += [(iatm, jatm)]

        for three in three_atoms:
            for pair in it.combinations(three, 2):
                iatm, jatm = min(pair), max(pair)
                pairs += [(iatm, jatm)]
                
        for four in four_atoms:
            for pair in it.combinations(four, 2):
                iatm, jatm = min(pair), max(pair)
                pairs += [(iatm, jatm)]

        for four in four_atoms_imp:
            for pair in it.combinations(four, 2):
                iatm, jatm = min(pair), max(pair)
                pairs += [(iatm, jatm)]

        return sorted(set(pairs), key=lambda x:x[0])

    def _make_bonded14_pairs(self):

        pairs = []
        four_atoms = self.get_torsion_info()['four_atoms']
        for four in four_atoms:
            pair = four[0], four[3]
            iatm, latm = min(pair), max(pair)
            pairs += [(iatm, latm)]

        three_pairs = []
        three_atoms = self.get_angle_info()['three_atoms']
        for three in three_atoms:
            pair = three[0], three[2]
            iatm, katm = min(pair), max(pair)
            three_pairs += [(iatm, katm)]

        return sorted(set(pairs)-set(three_pairs), key=lambda x:x[0])

    #TODO
    def _make_bonded14_pairs_old(self):

        cur_natom = 0
        two_atoms = []

        # loop for the number of the kinds for molcule
        for minfo in self.get_mol_info():

            tor_info = self.get_topology().get_torsion_info(minfo['name'])
            if tor_info:

                # loop for the number of molecules
                for imol in range( minfo['nmol'] ):

                    atoms_flags = zip(
                            tor_info['nonbonded_flags'],
                            tor_info['four_atoms'] )

                    two_atoms.extend(
                        [ (iatoms[0]+cur_natom, iatoms[3]+cur_natom) 
                            for flag, iatoms in atoms_flags if flag==1 ] )

                    cur_natom += minfo['natom']

            else:
                cur_natom += minfo['natom'] * minfo['nmol']

        self._bonded14_pairs = two_atoms

    def _get_idx_to_ipair(self, bonded_type):
        dictname = 'i'+bonded_type+'_to_ipair'
        to_ipair = getattr(self, '_'+dictname)
        if to_ipair is None:
            make_method = getattr(self, '_make_'+dictname)
            to_ipair = make_method()
        return to_ipair

    def _apply_target_to_bondtype(self, iatoms_list, info):

        target_atoms = self.__target_atoms

        new_infos = {key:[] for key in info.keys()}
        for i_1, iatoms in enumerate(iatoms_list):
            if not self._contains(target_atoms, iatoms): continue

            for key, new_info in new_infos.items():
                new_infos[key].append( info[key][i_1] )

        for key, new_info in new_infos.items():
            info[key] = new_info

    def _contains(self, target_atoms, iatoms):
        for iatm in iatoms:
            if iatm in target_atoms: return True
        else:
            return False

    def _make_ibnd_to_ipair(self):
        iatoms_list  = self.get_bond_info()['two_atoms']
        bonded_pairs = self.get_bonded_pairs()

        # make table
        ibnd_to_ipair = lib_bonded_pair.get_ibnd_to_ipair(
                numpy.array(iatoms_list), numpy.array(bonded_pairs) )

        return ibnd_to_ipair

    def _make_iang_to_ipair(self):
        iatoms_list  = self.get_angle_info()['three_atoms']
        bonded_pairs = self.get_bonded_pairs()

        # make table
        if len(iatoms_list) == 0:
            iang_to_ipair = numpy.array([[]])
        else:
            iang_to_ipair = lib_bonded_pair.get_iang_to_ipair(
                    numpy.array(iatoms_list), numpy.array(bonded_pairs) )

        return iang_to_ipair

    def _make_itor_to_ipair(self):
        iatoms_list  = self.get_torsion_info()['four_atoms']
        bonded_pairs = self.get_bonded_pairs()

        # make table
        if len(iatoms_list) == 0:
            itor_to_ipair = numpy.array([[]])
        else:
            itor_to_ipair = lib_bonded_pair.get_itor_to_ipair(
                    numpy.array(iatoms_list), numpy.array(bonded_pairs) )

        return itor_to_ipair

    def _make_iimp_to_ipair(self):
        iatoms_list  = self.get_improper_info()['four_atoms']
        bonded_pairs = self.get_bonded_pairs()

        # make table
        if len(iatoms_list) == 0:
            iimp_to_ipair = numpy.array([[]])
        else:
            iimp_to_ipair = lib_bonded_pair.get_itor_to_ipair(
                    numpy.array(iatoms_list), numpy.array(bonded_pairs) )

        return iimp_to_ipair

    def _make_i14_to_ipair(self):
        iatoms_list  = self.get_bonded14_pairs()
        bonded_pairs = self.get_bonded_pairs()

        # make table
        if len(iatoms_list) == 0:
            i14_to_ipair = numpy.array([[]])
        else:
            i14_to_ipair = lib_bonded_pair.get_ibnd_to_ipair(
                    numpy.array(iatoms_list), numpy.array(bonded_pairs) )

        return i14_to_ipair


class ConverterBase(ConverterPrintable, TableMaker):
    __metaclass__ = ABCMeta

    def __init__(self, topology):
        self.__topology = topology
        TableMaker.__init__(self)

    @abstractmethod
    def convert(self):
        """Convert program-specific data to the CURP data structure."""
        return

    @abstractmethod
    def get_natom(self):
        """Get the number of atoms."""
        return

    @abstractmethod
    def get_mol_info(self):
        """Get the molecular information."""
        return

    @abstractmethod
    def get_atom_info(self):
        """Get the dictionary of atom ids, names, masses, elems, vdw_radii."""
        return

    @abstractmethod
    def get_residue_info(self):
        """Get the dictionary of residue ids and names for each atom."""
        return

    @abstractmethod
    def get_bond_info(self):
        """Get bond dictionary of two_atoms, force_consts, and length_eqs."""
        return

    @abstractmethod
    def get_angle_info(self):
        """Get angle dictionary of three_atoms, force_consts, and theta_eqs."""
        return

    @abstractmethod
    def get_torsion_info(self):
        """Get torsion dictionary of four_atoms, num_torsions, num_freqs,
        force_consts, and initial_phases.
        """
        return

    @abstractmethod
    def get_improper_info(self):
        """Get improper torsion dictionary of four_atoms, num_torsions,
        num_freqs, force_consts, and initial_phases."""
        return

    @abstractmethod
    def get_coulomb_info(self):
        """Get coulomb dictionary of charges."""
        return

    @abstractmethod
    def get_vdw_info(self):
        """Get vdw dictionary of atom_types, c6s, c12s."""
        return

    @abstractmethod
    def get_pbc_info(self):
        """Get the box for beriodic boundary condition."""
        return

    def get_topology(self):
        return self.__topology


