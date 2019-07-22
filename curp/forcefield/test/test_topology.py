"""Test of curp forcefield topology.

To run: $ nosetests topology.py
"""
import parmed
from nose.tools import eq_
from curp.forcefield import topology


class TestIds:
    def setup(self):
        self.inter = topology.Interactions()
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

        # H-C-N-H CL
        # 1 0 2 3

        s.add_atom(parmed.Atom(atomic_number=6), 'ALA', 1, 'A')
        s.add_atom(parmed.Atom(atomic_number=1), 'ALA', 1, 'A')
        s.add_atom(parmed.Atom(atomic_number=7), 'ALA', 1, 'A')
        s.add_atom(parmed.Atom(atomic_number=1), 'ALA', 1, 'A')
        s.add_atom(parmed.Atom(atomic_number=17), 'CLA', 1, 'B')

        bond = list(range(3))
        angle = list(range(2))
        bond[0] = parmed.Bond(s.atoms[0], s.atoms[1], parmed.BondType(10., 1.))
        bond[1] = parmed.Bond(s.atoms[0], s.atoms[2], parmed.BondType(10., 1.))
        bond[2] = parmed.Bond(s.atoms[2], s.atoms[3], parmed.BondType(10., 1.))

        angle[0] = parmed.Angle(s.atoms[1], s.atoms[0], s.atoms[2],
                              parmed.AngleType(10., 166.))
        angle[1] = parmed.Angle(s.atoms[0], s.atoms[2], s.atoms[3],
                              parmed.AngleType(5., 66.))

        dihedral1 = parmed.Dihedral(s.atoms[0], s.atoms[1], s.atoms[2],
                                    s.atoms[3])
        dihedral1.type = parmed.DihedralType(phi_k=0.000, per=1, phase=0.000,
                                             scee=1.000, scnb=1.000)

        for b in bond: s.bonds.append(b)
        for a in angle: s.angles.append(a)
        s.dihedrals.append(dihedral1)
        self.tpl = topology.Topology(s)

    def test_create_topology(self):
        eq_(self.tpl.bonds.ids.tolist(), [[1, 2], [1, 3], [3,4]])
        eq_(self.tpl.angles.ids.tolist(), [[2, 1, 3], [1, 3, 4]])
        eq_(self.tpl.dihedrals.ids.tolist(), [[1, 2, 3, 4]])

    def test_forcefield_constants(self):
        eq_(self.tpl.bonds.k, [10., 10., 10.])
        eq_(self.tpl.bonds.req, [1., 1., 1.])
        eq_(self.tpl.angles.k, [10., 5.])
        eq_(self.tpl.angles.theteq, [166., 66.])

    def test_nonbonded_table(self):
        eq_(list(self.tpl.nonbonded_table), [[1, 5, 5],
                                             [2, 5, 5],
                                             [3, 5, 5],
                                             [4, 5, 5]])
