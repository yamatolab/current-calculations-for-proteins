"""
Topology and Interactions objects.


To understand @property and @x.setter:
https://docs.python.org/3/library/functions.html?highlight=property#property
"""
from __future__ import print_function

import os
import sys
import itertools

import numpy
from parmed.topologyobjects import (Bond, Angle, Dihedral, Improper,
                                    UreyBradley, Cmap, TrackedList)

# curp package
from curp.table.interact_table import InteractionTable
import curp.clog as logger
from curp.forcefield import lib_bonded_pair


class Interactions(object):
    """Interactions representation

    Parameters
    ----------
    parmed : parmed.topologyobjects.TrackedList
        TrackedList containing parmed.topologyobjects.Bond, .Angle,
        .Dihedral, .Improper, .UreyBradley or .Cmap
    ids: list of list of int
        Lists of ids of the atoms per interactions for fortran.
        (atom ids starting by 1, not 0).

    Attributes
    ----------
    parmed : parmed.topologyobjects.TrackedList
             See Parameters description.
    ids: list of list of int
         See Parameters description.
         If parmed used, all atoms ids are shifted by 1
    num_atm : int
        Number of atoms in one interaction
        Ex: 3 for angles or urey-bradleys
    pairs : list of tuples of int
        List of pair of atoms interacting.
    to_ipair : list of list of int or None
        For each interaction index, returns the associated pairs index.
        Created and edited only when pairs are edited.
        For fortran (pairs indexes starting by 1).
    ff_cst : dict
        Each key is an attribute relevant for forcefield calculations
        automatically extracted from parmed[0].type.
        Each key is also a class attribute.

    Methods
    -------
    make_to_ipair(pairs)
        Set the to_ipair attribute depending on ids attribute and pairs.
    flog()
        Called by clog.info().

    Examples
    --------
    >>> from curp.forcefield.topology import Interactions
    >>> angles = Interactions()
    >>> angles.ids = [[1, 2, 3],
    ...               [4, 1, 3]]
    >>> print(angles.pairs)
    [(1, 2), (1, 3), (1, 4), (2, 3), (3, 4)]
    >>> angles.make_to_ipair(angles.pairs)
    >>> print(angles.to_ipair)
    array([[0, 1, 3],
           [-2, -4, 1]])
    """
    def __init__(self, parmed=[], ids=[], ff_cst={}):
        self._ids = []
        self.ff_cst = {}

        # If a parmed interactions TrackedList is given, ids and ff_cst
        # are automatically made.
        self.parmed = parmed

        # If a parmed TrackedList was entered, setting ids or ff_cst
        # will return an error.
        if ids:
            self.ids = ids
        if ff_cst:
            self.ff_cst = ff_cst

        self._to_ipair = None

    def __str__(self):
        """
        Representation of the class using print()
        """
        ff_cst_names = list(self.ff_cst.keys())
        ff_list = [self.ff_cst[key] for key in ff_cst_names]

        # Create the format
        fmt = '{:>5} {:>20}'
        for key in ff_cst_names: fmt += ' {:>12}'
        fmt += '\n'

        # Create the header
        rpr = fmt.format('id', 'iatoms', *ff_cst_names)

        # Create the body
        fmt = fmt.replace(':>12', ':>12.7f')
        for inter_id, infos in enumerate(zip(self.ids.tolist(),
                                             *ff_list), 1):
            rpr += fmt.format(inter_id, str(infos[0]), *infos[1:])
        return rpr

    def __repr__(self):
        """
        Representation of the class using repr()
        """
        rpr = ('{!s}(\n'
               '    ids={!r},\n'
               '    ff_cst={!r})').format(self.__class__,
                                          self.ids,
                                          self.ff_cst)
        return rpr

    # To understand @property see property decorator.
    @property
    def parmed(self):
        return self._parmed

    @parmed.setter
    def parmed(self, new_parmed):
        self._parmed = new_parmed
        self._make_ff_cst_lists()
        self._num_atm = self._make_num_atm()
        self._ids = self._make_ids()
        self._pairs = self._make_pairs()

    @property
    def ids(self):
        return self._ids

    @ids.setter
    def ids(self, new_ids):
        if self.parmed == []:
            self._ids = new_ids
            self._num_atm = self._make_num_atm()
        else:
            raise Exception('ids can not be edited if parmed attribute exists')

        self._pairs = self._make_pairs()

    @property
    def num_atm(self):
        return self._num_atm

    @property
    def pairs(self):
        return self._pairs

    @property
    def to_ipair(self):
        return self._to_ipair

    def _make_ff_cst_lists(self):
        """Make forcefield constant lists and set them as attributes"""
        if len(self.parmed) > 0:
            for key in self.parmed[0].type.__dict__:
                if key not in ('used', 'penalty', 'idx', 'list', '_idx'):
                    ff_cst_list = [getattr(inter.type, key)
                                   for inter in self.parmed]
                    setattr(self, key, ff_cst_list)
                    self.ff_cst[key] = getattr(self, key)

    def _make_num_atm(self):
        """Return the number of atoms in one interaction"""
        num_atm = 1
        try:
            while 'atom{}'.format(num_atm+1) in dir(self.parmed[0]):
                num_atm += 1
        except IndexError:
            try:
                num_atm = len(self.ids[0])
            except AttributeError:
                num_atm = 0
            except IndexError:
                num_atm = 0
        return(num_atm)

    def _make_ids(self):
        """Return lists of atom ids for each interaction"""
        ids = numpy.zeros((len(self.parmed), self.num_atm), int)
        for i, interaction in enumerate(self.parmed):
            for j in range(self.num_atm):
                atom = getattr(interaction, 'atom{}'.format(j+1))
                ids[i, j] = atom.idx
        ids = ids + 1
        return(ids)

    def _make_pairs(self):
        """Return list of all sorted unique interacting atom pairs"""
        pairs = []   # ipair => (iatm, jatm)
        for atom_ids in self.ids:
            for pair in itertools.combinations(atom_ids, 2):
                iatm, jatm = min(pair), max(pair)
                pairs.append((iatm, jatm))
        return sorted(set(pairs), key=lambda x: x[0])

    def make_to_ipair(self, pairs):
        """Link interaction table with bonded_pairs table.

        Set self.to_ipair given a list of atom id pairs and self.ids.

        Parameters
        ----------
        pairs: list of list of integers
        """

        # Compute the number of combinations between atoms in one interaction
        from math import factorial
        num_comb = (factorial(self.num_atm)
                    // (factorial(2)*factorial(self.num_atm - 2)))
        to_ipair = numpy.zeros([len(self.ids), num_comb], int)

        # Find pair id for each pair combinations of an interaction
        for inter, iatoms in enumerate(self.ids):
            for icomb, pair in enumerate(itertools.combinations(iatoms, 2)):
                sorted_pair = (min(pair), max(pair))
                ipair = pairs.index(sorted_pair)
                if sorted_pair != pair:
                    ipair = - ipair
                to_ipair[inter, icomb] = ipair

        self._to_ipair = to_ipair


class TableMaker(object):
    """Define required methods for making interaction table"""
    def __init__(self):
        self._bonded_inter = {}   # Bonded Interactions objects
        self._nonbonded_inter = {}    # Non-bonded Interactions objects
        self._bonded_pairs = None
        self._nonbonded_table = None
        self._natom = 0

    @property
    def bonded_inter(self):
        return self._bonded_inter

    @property
    def bonded_pairs(self):
        if self._bonded_pairs is None:
            self._bonded_pairs = self._concat_pairs(self.bonded_inter)
        return self._bonded_pairs

    @property
    def nonbonded_table(self):
        if self._nonbonded_table is None:
            self._nonbonded_table = self._make_nonbonded_table()
        return self._nonbonded_table

    @property
    def natom(self):
        return self._natom

    def _make_nonbonded_table(self):
        """Create non-bonded InteractionTable

        Returns
        -------
        new_table : curp.table.interact_table.InteractionTable
            form: [[iatm, jatm_beg, jatm_end]]
            with [jatm_beg, jatm_end] a range of atoms with which
            iatm is not interacting.

        Examples
        --------
        Considering a topology file containing two water molecules:
        H-O-H H-O-H
        1 2 3 4 5 6

        list(new_table) would be:
        [[1, 4, 6],
         [2, 4, 6],
         [3, 4, 6]]
        """
        bonded_pairs = numpy.array(self.bonded_pairs)
        int_table = InteractionTable(self.natom)

        from curp.forcefield import lib_nonbond
        mod_tab = lib_nonbond.without_bonded
        mod_tab.setup(bonded_pairs, self.natom)

        int_table = numpy.array(list(int_table))
        nonbonded_table, ntable = mod_tab.get_nonbonded_table(int_table)
        new_table = InteractionTable(
                base_table=nonbonded_table[:ntable].tolist())

        return new_table

    def _concat_pairs(self, inter_lists):
        """Concatenate pairs from a dictionary of Interactions instances.

        Parameters
        ----------
        inter_lists : dict of Interactions

        Returns
        -------
        list_pairs : list of list of ints
            Sorted, unique pairs of atom ids.
        """

        lists_pairs = []
        for key, inter in inter_lists.items():
            lists_pairs += inter.pairs
        return sorted(set(lists_pairs))


class Topology(TableMaker):
    def __init__(self, parmed):
        TableMaker.__init__(self)

        self.parmed = parmed
        self._natom = len(self.parmed.atoms)
        self.bonds = Interactions(self.parmed.bonds)
        self.angles = Interactions(self.parmed.angles)
        self.dihedrals = Interactions(self.parmed.dihedrals)
        self.impropers = Interactions(self.parmed.impropers)

        self._bonded_inter = {
                'bonds': self.bonds,
                'angles': self.angles,
                'dihedrals': self.dihedrals,
                'impropers': self.impropers
                }

    def get_decomp_list(self, dtype='all'):
        """Get the list that decomposes all potential.
        dtype = 'all', 'bonded', 'bonded14', 'bonded+' (bonded and 14),
        'nonbonded' or 'nonbonded+' (nonbonded and 14)
        """
        bonded_list = ['bond', 'angle', 'torsion', 'improper']
        bonded14_list = ['coulomb14', 'vdw14']
        nonbonded_list = ['coulomb', 'vdw']

        if dtype == 'bonded': return bonded_list
        elif dtype == 'bonded14': return bonded14_list
        elif dtype == 'bonded+': return bonded_list + bonded14_list
        elif dtype == 'nonbonded': return nonbonded_list
        elif dtype == 'nonbonded+': return nonbonded_list + bonded14_list
        else: return bonded_list + bonded14_list + nonbonded_list

    # For amber only:

    def _make_bonded14_pairs(self):
        for four in four_atoms:
            pair = four[0], four[3]
            iatm, latm = min(pair), max(pair)
            pairs += [(iatm, latm)]

        three_pairs = []
        three_atoms = [[angle.atom1.idx, angle.atom2.idx, angle.atom3.idx]
                       for angle in self.tpl.angles]
        for three in three_atoms:
            pair = three[0], three[2]
            iatm, katm = min(pair), max(pair)
            three_pairs += [(iatm, katm)]
        return sorted(set(pairs)-set(three_pairs), key=lambda x: x[0])

    @property
    def bonded14_pairs(self):
        if self._bonded14_pairs is None:
            self._bonded14_pairs = self._concat_pairs(self.bonded14_inter)
        return self._bonded14_pairs
