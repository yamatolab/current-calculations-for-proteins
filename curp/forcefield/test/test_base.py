"""Test of curp forcefield base.

To run: $ nosetests base.py
"""
import parmed
from nose.tools import eq_
from curp.forcefield import base


class TestIds:
    def setup(self):
        self.inter = base.Interactions()
        self.inter.ids = [[1, 2, 3],
                          [4, 1, 3]]

    def test_num_atm(self):
        eq_(self.inter.num_atm, 3)

    def test_pairs(self):
        eq_(self.inter.pairs,
            [(1, 2), (1, 3), (1, 4), (2, 3), (3, 4)])

    def test_make_ids_to_pairs(self):
        self.inter.make_to_ipair(self.inter.pairs)
        eq_(self.inter.to_ipair.tolist(),
            [[0, 1, 3], [-2, -4, 1]])


class TestParmedTopology:
    def setup(self):
        s = self.s = parmed.structure.Structure()
        s.add_atom(parmed.Atom(atomic_number=6), 'ALA', 1, 'A')
        s.add_atom(parmed.Atom(atomic_number=1), 'ALA', 1, 'A')
        s.add_atom(parmed.Atom(atomic_number=1), 'ALA', 1, 'A')
        bond1 = parmed.Bond(s.atoms[0], s.atoms[1], parmed.BondType(10., 1.))
        bond2 = parmed.Bond(s.atoms[0], s.atoms[2], parmed.BondType(10., 1.))
        angle1 = parmed.Angle(s.atoms[1], s.atoms[0], s.atoms[2],
                              parmed.AngleType(10., 180.))
        s.bonds.append(bond1)
        s.bonds.append(bond2)
        s.angles.append(angle1)
        self.tpl = base.Topology(s)

    def test_create_topology(self):
        eq_(self.tpl.bonds.ids.tolist(), [[0, 1], [0, 2]])
        eq_(self.tpl.angles.ids.tolist(), [[1, 0, 2]])

    def test_nonbonded_table(self):
        eq_(self.tpl.nonbonded_table, 1)
