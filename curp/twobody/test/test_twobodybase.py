from nose.tools import *

import os, sys
import numpy
import parmed

from curp.forcefield import Topology
from curp.twobody import TwoBodyForce

class DummySetting:
    class Curp:
        coulomb_cutoff_length = 10.0
        vdw_cutoff_length = 10.0
    curp = Curp()

def get_topology():
    """
    Create a parmed Structure such as:
    H-C-N-H CL
    1 0 2 3 4
    """
    s = parmed.structure.Structure()

    # Create atoms
    s.add_atom(parmed.Atom(name='C', atomic_number=6), 'ALA', 1, 'A')
    s.add_atom(parmed.Atom(name='HC', atomic_number=1), 'ALA', 1, 'A')
    s.add_atom(parmed.Atom(name='N', atomic_number=7), 'ALA', 1, 'A')
    s.add_atom(parmed.Atom(name='HN', atomic_number=1), 'ALA', 1, 'A')
    s.add_atom(parmed.Atom(name='CL', atomic_number=17), 'CLA', 1, 'B')

    # Create interactions
    bond = list(range(3))
    angle = list(range(2))
    bond[0] = parmed.Bond(s.atoms[0], s.atoms[1], parmed.BondType(k=10., req=1.))
    bond[1] = parmed.Bond(s.atoms[0], s.atoms[2], parmed.BondType(10., 1.))
    bond[2] = parmed.Bond(s.atoms[2], s.atoms[3], parmed.BondType(10., 1.))

    angle[0] = parmed.Angle(s.atoms[1], s.atoms[0], s.atoms[2],
                            parmed.AngleType(k=10., theteq=166.))
    angle[1] = parmed.Angle(s.atoms[0], s.atoms[2], s.atoms[3],
                            parmed.AngleType(5., 66.))

    dihedral1 = parmed.Dihedral(s.atoms[0], s.atoms[1], s.atoms[2],
                                s.atoms[3])
    dihedral1.type = parmed.DihedralType(phi_k=0.000, per=1, phase=0.000,
                                            scee=1.000, scnb=1.000)

    for b in bond: s.bonds.append(b)
    for a in angle: s.angles.append(a)
    s.dihedrals.append(dihedral1)

    tpl = Topology(s)
    return(tpl)

class TestTwoBodyForce:
    def setup(self):
        """Get TwoBodyForce from its Topology object."""
        tpl = get_topology()
        setting = DummySetting()
        self.tbcal = TwoBodyForce(tpl, setting)
        self.interaction_table = []

    def test_setup(self):
        self.tbcal.setup(self.interaction_table)


# def cal_bond(r_i, r_j ,l_eq, kb):
#     l2 = numpy.dot(r1, r2)
#     return kb * l2
#
# class TestAmber99:
#
#     def test_bond(self):
#
#         two_atoms    = [(1,2), (1,3), (4,5)]
#         # initialize with coordinate
#         crd = numpy.array(
#             [[1.0, 2.0, 3.0],
#              [2.0, 1.0, 1.5],
#              [3.0, 3.0, 3.0],
#              [0.0, 0.0, 0.0],
#              [5.5, 5.1, 4.5]])
#         tbcal.initialize(crd)
#
#         res = tbcal.cal_bonded('bond')
#
#         zero_tbfs = res['zero_tbforces']
#         bonded_tbfs = res['tbforces']
#
#         eq_(self.get_nonzero_pairs(zero_tbfs), set(two_atoms))
#
#     def get_bond_info(self):
#         return dict(
#             two_atoms    = [(1,2), (1,3), (4,5)],
#             force_consts = [1.0, 1.5, 1.9],
#             length_eqs   = [1.0, 1.2, 0.9] )
#
#         # right values
#         info = tpl.get_bond_info()
#         a_res = []
#         for pair, kb, l_eq in zip(info['two_atoms'],
#                                   info['force_consts'],
#                                   info['length_eqs']):
#             iatm, jatm = pair
#             r_i, r_j = crd[iatm-1], crd[jatm-1]
#             a_res.append( pair, cal_bond(r_i, r_j, l_eq, kb))
#
#         import numpy.testing as npt
#         n = len(a_res)
#         for (p, a_tbf), (iatm, jatm, tbf) in zip(a_res,
#                 self.gen_bonded_tbfs(bonded_tbfs, zero_tbfs)):
#
#             eq_((iatm, jatm), p)
#             npt.assert_almost_equal(tbf, a_tbf)
#
#     def get_nonzero_pairs(self, zero_tbfs):
#         natom = len(zero_tbfs)
#         n_ext = len(zero_tbfs[0])
#
#         bonded_pairs = []
#         for iatm_1 in range(natom):
#             iatm = iatm_1 + 1
#             for iext_1 in range(n_ext):
#                 iext = iext_1 + 1
#                 jatm = iatm + iext
#                 if zero_tbfs[iatm_1, iext_1]: continue
#                 bonded_pairs.append( (iatm, jatm) )
#
#         return set(bonded_pairs)
#
#     def gen_bonded_tbfs(self, tbforces, zero_tbfs):
#         natom = len(zero_tbfs)
#         n_ext = len(zero_tbfs[0])
#
#         for iatm_1 in range(natom):
#             iatm = iatm_1 + 1
#             for iext_1 in range(n_ext):
#                 iext = iext_1 + 1
#                 jatm = iatm + iext
#                 if zero_tbfs[iatm_1, iext_1]: continue
#                 yield iatm, jatm, tbforces[iatm_1,iext_1]
#
#
#
#         # for iatm_1, zero, tbf in enumerate(zip(zero_tbfs, bonded_tbfs)):
#         #     iatm = iatm_1 + 1
#         #     print(iatm,
#         # #     if iszero: continue
#         #     # print(tbf)
#         # # print(res['zero_tbforces'])
#         # eq_(iszeros[1,1], False)
#
#         # do iatm=1, natom
#         # do iext=1, nextent
#         #     jatm = iatm + iext
#         #     tbf_abs = dot_product(tbforces(iatm,iext,:), tbforces(iatm,iext,:))
#         #     if (-0.0001 < tbf_abs .and. tbf_abs < 0.0001) cycle
#         #     write(*, '(2I5,3F15.7)'), iatm, jatm, tbforces(iatm, iext, :)
#         # end do
#         # end do
#
#         # eq_(res, 1)
#
#     def test_angle(self):
#
#         # initialize with coordinate
#         crd = numpy.array(
#             [[1.0, 2.0, 3.0],
#              [2.0, 1.0, 1.5],
#              [3.0, 3.0, 3.0],
#              [0.0, 0.0, 0.0],
#              [5.5, 5.1, 4.5]])
#         tbcal.initialize(crd)
#
#         res = tbcal.cal_bonded('angle')
#
#     def test_torsion(self):
#
#         # initialize with coordinate
#         crd = numpy.array(
#             [[1.0, 2.0, 3.0],
#              [2.0, 1.0, 1.5],
#              [3.0, 3.0, 3.0],
#              [0.0, 0.0, 0.0],
#              [5.5, 5.1, 4.5]])
#         tbcal.initialize(crd)
#
#         res = tbcal.cal_bonded('torsion')
#
#     def test_improper(self):
#
#         # initialize with coordinate
#         crd = numpy.array(
#             [[1.0, 2.0, 3.0],
#              [2.0, 1.0, 1.5],
#              [3.0, 3.0, 3.0],
#              [0.0, 0.0, 0.0],
#              [5.5, 5.1, 4.5]])
#         tbcal.initialize(crd)
#
#         res = tbcal.cal_bonded('improper')
#
#     def test_coulomb(self):
#         table = [[1, 2, 5], [2, 3, 5], [3,4,5], [4, 5, 5] ]
#         res = tbcal.cal_coulomb(table)
#         tbcal.output_nonbonded(res, 'coulomb', table)
#
#     def test_vdw(self):
#         table = [[1, 2, 5], [2, 3, 5], [3,4,5], [4, 5, 5] ]
#         res = tbcal.cal_vdw(table)
#         tbcal.output_nonbonded(res, 'vdw', table)

# # crd = numpy.array(
# #     [[2.0, 3.0, 4.1],
# #      [3.0, 2.0, 2.5],
# #      [4.0, 4.0, 4.0],
# #      [1.0, 1.0, 1.0],
# #      [9.0, 9.1, 9.2]])
# # tbcal.initialize(crd, check=True)
# # tbcal.setup()
# # results = tbcal.cal_bonded()
# # output(results)

# # from benchmarker import Benchmarker
# # with Benchmarker(width=20) as bm:

# #     with bm('setup'):
# #         crd = numpy.array(
# #             [[2.0, 3.0, 4.1],
# #              [3.0, 2.0, 2.5],
# #              [4.0, 4.0, 4.0],
# #              [1.0, 1.0, 1.0],
# #              [9.0, 9.1, 9.2]])
# #         tbcal.initialize(crd, check=False)
# #         tbcal.setup()

# #     with bm('bonded'):
# #         results = tbcal.cal_bonded()
# #         output(results)

# #     table = [[1, 2, 5], [2, 3, 3], [2, 5, 5]]
# #     with bm('coulomb'):
# #         results = tbcal.cal_coulomb(table)
# #         output(results)

# #     with bm('vdw'):
# #         results = tbcal.cal_vdw(table)
# #         output(results)


