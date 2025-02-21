from nose.tools import *


import os, sys
import numpy as np
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)


from forcefield.amber99 import ConverterBase
from twobody.amber99 import TwoBodyForce

################################################################################
class DummySetting:
    class Curp:
        coulomb_cutoff_length = 10.0
        vdw_cutoff_length = 10.0
    curp = Curp()

################################################################################
class DummyTopology(ConverterBase):

    def __init__(self):
        ConverterBase.__init__(self, None)

    def get_pbc_info(self): pass
    def get_mol_info(self): pass
    def convert(self): pass
    def get_residue_info(self): pass

    def get_natom(self):
        return 5

    def get_bond_info(self):
        return dict(
            two_atoms    = [(1,2), (1,3), (4,5)],
            force_consts = [1.0, 1.5, 1.9],
            length_eqs   = [1.0, 1.2, 0.9] )

    def get_angle_info(self):
        return dict(
            three_atoms  = [(1,2,3), (2,3,4), (1,2,5)],
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


tpl = DummyTopology()
tbcal = TwoBodyForce(tpl, DummySetting())
tbcal.setup()

def cal_bond(r_i, r_j ,l_eq, kb):
    l2 = np.dot(r1, r2)
    return kb * l2

# def cal_angle(r_i, r_j, r_k, theta_eq, kt):





################################################################################
class TestAmber99:

    def test_bond(self):

        two_atoms    = [(1,2), (1,3), (4,5)]
        # initialize with coordinate
        crd = np.array(
            [[1.0, 2.0, 3.0],
             [2.0, 1.0, 1.5],
             [3.0, 3.0, 3.0],
             [0.0, 0.0, 0.0],
             [5.5, 5.1, 4.5]])
        tbcal.initialize(crd)

        res = tbcal.cal_bonded('bond')
        
        zero_tbfs = res['zero_tbforces']
        bonded_tbfs = res['tbforces']

        eq_(self.get_nonzero_pairs(zero_tbfs), set(two_atoms))

    def get_bond_info(self):
        return dict(
            two_atoms    = [(1,2), (1,3), (4,5)],
            force_consts = [1.0, 1.5, 1.9],
            length_eqs   = [1.0, 1.2, 0.9] )

        # right values
        info = tpl.get_bond_info()
        a_res = []
        for pair, kb, l_eq in zip(info['two_atoms'],
                                  info['force_consts'],
                                  info['length_eqs']):
            iatm, jatm = pair
            r_i, r_j = crd[iatm-1], crd[jatm-1]
            a_res.append( pair, cal_bond(r_i, r_j, l_eq, kb))

        import numpy.testing as npt
        n = len(a_res)
        for (p, a_tbf), (iatm, jatm, tbf) in zip(a_res,
                self.gen_bonded_tbfs(bonded_tbfs, zero_tbfs)):

            eq_((iatm, jatm), p)
            npt.assert_almost_equal(tbf, a_tbf)

    def get_nonzero_pairs(self, zero_tbfs):
        natom = len(zero_tbfs)
        n_ext = len(zero_tbfs[0])

        bonded_pairs = []
        for iatm_1 in range(natom):
            iatm = iatm_1 + 1
            for iext_1 in range(n_ext):
                iext = iext_1 + 1
                jatm = iatm + iext
                if zero_tbfs[iatm_1, iext_1]: continue
                bonded_pairs.append( (iatm, jatm) )

        return set(bonded_pairs)

    def gen_bonded_tbfs(self, tbforces, zero_tbfs):
        natom = len(zero_tbfs)
        n_ext = len(zero_tbfs[0])

        for iatm_1 in range(natom):
            iatm = iatm_1 + 1
            for iext_1 in range(n_ext):
                iext = iext_1 + 1
                jatm = iatm + iext
                if zero_tbfs[iatm_1, iext_1]: continue
                yield iatm, jatm, tbforces[iatm_1,iext_1]



        # for iatm_1, zero, tbf in enumerate(zip(zero_tbfs, bonded_tbfs)):
        #     iatm = iatm_1 + 1
        #     print(iatm, 
        # #     if iszero: continue
        #     # print(tbf)
        # # print(res['zero_tbforces'])
        # eq_(iszeros[1,1], False)

        # do iatm=1, natom
        # do iext=1, nextent
        #     jatm = iatm + iext
        #     tbf_abs = dot_product(tbforces(iatm,iext,:), tbforces(iatm,iext,:))
        #     if (-0.0001 < tbf_abs .and. tbf_abs < 0.0001) cycle
        #     write(*, '(2I5,3F15.7)'), iatm, jatm, tbforces(iatm, iext, :)
        # end do
        # end do

        # eq_(res, 1)

    def test_angle(self):

        # initialize with coordinate
        crd = np.array(
            [[1.0, 2.0, 3.0],
             [2.0, 1.0, 1.5],
             [3.0, 3.0, 3.0],
             [0.0, 0.0, 0.0],
             [5.5, 5.1, 4.5]])
        tbcal.initialize(crd)

        res = tbcal.cal_bonded('angle')

    def test_torsion(self):

        # initialize with coordinate
        crd = np.array(
            [[1.0, 2.0, 3.0],
             [2.0, 1.0, 1.5],
             [3.0, 3.0, 3.0],
             [0.0, 0.0, 0.0],
             [5.5, 5.1, 4.5]])
        tbcal.initialize(crd)

        res = tbcal.cal_bonded('torsion')

    def test_improper(self):

        # initialize with coordinate
        crd = np.array(
            [[1.0, 2.0, 3.0],
             [2.0, 1.0, 1.5],
             [3.0, 3.0, 3.0],
             [0.0, 0.0, 0.0],
             [5.5, 5.1, 4.5]])
        tbcal.initialize(crd)

        res = tbcal.cal_bonded('improper')

    def test_coulomb(self):
        table = [[1, 2, 5], [2, 3, 5], [3,4,5], [4, 5, 5] ]
        res = tbcal.cal_coulomb(table)
        tbcal.output_nonbonded(res, 'coulomb', table)

    def test_vdw(self):
        table = [[1, 2, 5], [2, 3, 5], [3,4,5], [4, 5, 5] ]
        res = tbcal.cal_vdw(table)
        tbcal.output_nonbonded(res, 'vdw', table)

# # crd = np.array(
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
# #         crd = np.array(
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


