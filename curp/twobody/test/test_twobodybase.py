"""Test of curp basic TwoBodyForce object.

Requires Topology object to pass tests.

Don't forget to reinstall curp ($ pip install . in main directory,
where setup.py can be found) if any code modifications have been
made.
To run:
$ nosetests test_twobodybase.py
To print standard outputs anyways:
$ nosetests --nocapture test_twobodybase.py
or
$ nosetests -s test_twobodybase.py
"""
import os, sys
import numpy as np
from numpy.testing import assert_array_almost_equal
import parmed

from curp.forcefield import Topology
from curp.twobody import TwoBodyForce


class DummySetting:
    class Curp:
        coulomb_cutoff_length = 10.0
        vdw_cutoff_length = 10.0
    curp = Curp()


class TestTwoBodyForce:
    @classmethod
    def setup_class(cls):
        """Get TwoBodyForce from its Topology object."""
        cls.tpl = cls.make_topology()
        setting = DummySetting()
        cls.tbcal = TwoBodyForce(cls.tpl, setting)
        cls.interaction_table = []
        cls.tbcal.setup([])
        crd = np.array(
             [[0., 0., 0.],
              [0., 1., 0.],
              [2., 1., 0.],
              [2.+np.sqrt(3), 3., 3.],
              [5.5, 5.1, 4.5]])

        cls.test_tbforces = {}
        cls.test_tbforces['bond'] = np.array(
                                   [[0., 0., 0.],
                                    [0., 0., 0.],
                                    [0., 0., 0.],
                                    [2., 0., 0.],
                                    [0., 0., 0.],
                                    [1.5*np.sqrt(3), 3., 4.5]])

        cls.test_tbforces['angle'] = np.array(
                                   [[0., 0., 0.],
                                    [0., 0., 0.],
                                    [0., 0., 0.],
                                    [-0.46354444, 0., 0.],
                                    [0.46354444, 0.24841272, 0.37261908],
                                    [-0.26170911, -0.30219565, -0.45329347]])

        cls.test_tbforces['dihedral'] = np.array(
                                   [[0., 0.31450797, 0.],
                                    [1.03508466, 0.51754233, 0.],
                                    [-1.03508466, -0.83205029, -0.83205029],
                                    [-0.89717549, 0., 0.],
                                    [0.89717549, 0.48079489, 0.72119234],
                                    [0.13790917, 0.15924379, 0.23886569]])
        cls.tbcal.module.initialize(crd)

    @classmethod
    def make_topology(cls):
        """Create a parmed Structure.

        H-C-N-H CL
        0 1 2 3 4
        """
        s = parmed.structure.Structure()

        # Create atoms
        s.add_atom(parmed.Atom(name='HC', atomic_number=1), 'ALA', 1, 'A')
        s.add_atom(parmed.Atom(name='C', atomic_number=6), 'ALA', 1, 'A')
        s.add_atom(parmed.Atom(name='N', atomic_number=7), 'ALA', 1, 'A')
        s.add_atom(parmed.Atom(name='HN', atomic_number=1), 'ALA', 1, 'A')
        s.add_atom(parmed.Atom(name='CL', atomic_number=17), 'CLA', 1, 'B')

        # Create interactions
        # Three bonds, two angles, one dihedral
        bond = list(range(3))
        angle = list(range(2))
        bond[0] = parmed.Bond(s.atoms[0], s.atoms[1], parmed.BondType(k=1., req=1.))
        bond[1] = parmed.Bond(s.atoms[1], s.atoms[2], parmed.BondType(1., 1.))
        bond[2] = parmed.Bond(s.atoms[2], s.atoms[3], parmed.BondType(1., 1.))

        angle[0] = parmed.Angle(s.atoms[0], s.atoms[1], s.atoms[2],
                                parmed.AngleType(k=1., theteq=90.))
        angle[1] = parmed.Angle(s.atoms[1], s.atoms[2], s.atoms[3],
                                parmed.AngleType(1., 90.))

        dihedral1 = parmed.Dihedral(s.atoms[0], s.atoms[1], s.atoms[2],
                                    s.atoms[3])
        dihedral1.type = parmed.DihedralType(phi_k=1., per=1, phase=0.,
                                                scee=1., scnb=1.)

        for b in bond: s.bonds.append(b)
        for a in angle: s.angles.append(a)
        s.dihedrals.append(dihedral1)

        tpl = Topology(s)
        return(tpl)

    def tbforces_to_forces(self, tbforces, shape):
        """Sum twobody forces to obtain forces applied on each atom.

        Parameters
        ----------
        tbforces : 2D numpy.array
            Twobody Forces array as outputed by TwoBodyForces.cal_forces().
        shape : list of int
            Shape of the returned array.

        Returns
        -------
        forces : 2D numpy.array
            Forces deduced from summed twobody forces.

        Examples
        --------
        Only two atoms -> only one pair of atoms (1, 2).
        >>> tbforces = np.array([[1, 2, 3]])
        >>> tbforces_to_forces(self, tbforces, [2, 3])
        np.array([[1, 2, 3],
            [-1, -2, -3]])

        Indeed forces applied by 2 on 1 are the opposite of the ones
        applied by 1 on 2 (Newton's third law of motion).
        """
        forces = np.zeros(shape)
        pair_ids=[]
        for i, pair in enumerate(self.tpl.bonded_pairs):
            # Substract 1 to pass from fortran to python numbering
            forces[pair[0] - 1] += tbforces[i]
            forces[pair[1] - 1] -= tbforces[i]
        return(forces)

    def assert_bonded(self, btype):
        """Assert if twobody forces of bonded terms are correct.

        Assert both if the sum of twobody forces is equal to the
        computed forces, and if the twobody forces haven't changed
        (as defined in class setup).
        """
        res = self.tbcal._cal_bonded(btype)

        forces = res['forces']
        tbforces = res['tbforces']
        interactions = self.tpl.bonded_inter[btype]
        test_forces = self.tbforces_to_forces(tbforces, forces.shape)

        print(btype, '\nTwobody Forces:\n', tbforces,
                '\nForces applied on atoms:\n', forces,
              '\nForces as sum of tb forces:\n', test_forces)
        assert_array_almost_equal(tbforces, self.test_tbforces[btype],
                                 err_msg=('twobody forces of {} term were '
                                         'changed.').format(btype))
        assert_array_almost_equal(forces, test_forces,
                                  err_msg=('Sum of {} twobody forces are '
                                           'different from forces.').format(btype))

    def test_bond(self):
        self.assert_bonded('bond')

    def test_angle(self):
        self.assert_bonded('angle')

    def test_dihedral(self):
        self.assert_bonded('dihedral')
