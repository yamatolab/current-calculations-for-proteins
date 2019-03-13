import os, sys
import numpy
from nose.tools import *
import numpy.testing as npt


topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)

def tensor(x, y):
    """Calculate 3 x 3 tensor."""
    return numpy.array([[x[0]*y[0], x[0]*y[1], x[0]*y[2]],
                        [x[1]*y[0], x[1]*y[1], x[1]*y[2]],
                        [x[2]*y[0], x[2]*y[1], x[2]*y[2]]])
    
def get_iatm_to_itars(target_atoms, natom):
    iatm_to_itars = numpy.zeros( [natom], numpy.int)
    for itar_1, iatm in enumerate(target_atoms):
        iatm_to_itars[iatm-1] = itar_1 + 1
    return iatm_to_itars

def gen_kinetic_current(vel, masses, volumes, iatm_to_itars):
    for iatm_1, (v, m) in enumerate( zip(vel, masses) ):
        itar = iatm_to_itars[iatm_1]
        if itar != 0:
            vol  = volumes[itar-1]
            yield m * tensor(v, v) / vol

def gen_pair_to_tbfs(zero_tbfs, tbfs):
    natom = len(zero_tbfs)
    n_ext = len(zero_tbfs[0])

    for iatm_1 in range(natom):
        iatm = iatm_1 + 1
        for i_ext_1 in range(n_ext):
            i_ext = i_ext_1 + 1
            if zero_tbfs[iatm_1, i_ext_1]: continue
            jatm = iatm + i_ext
            yield iatm, jatm, tbfs[iatm_1, i_ext_1]

# Need for bonded current calculation
def get_extended_tbfs(pair_to_tbfs, natom, n_ext):
    zero_tbfs = numpy.ones((natom, n_ext),  numpy.bool)
    tbforces  = numpy.zeros((natom, n_ext, 3), numpy.float)
    for iatm, jatm, tbf in pair_to_tbfs:
        iatm_1 = iatm-1
        iext_1 = jatm - iatm - 1
        zero_tbfs[iatm_1, iext_1] = False
        tbforces[iatm_1, iext_1] = tbf

    return tbforces, zero_tbfs

# Need for nonbonded current calculation
def old_gen_table_tbfs(pair_to_tbfs, table, natom):

    # pair_to_tbfs = [
    #         (1, 2, numpy.array([0.2, 0.3, 0.4])),
    #         (1, 3, numpy.array([0.5, 0.6, 0.9])),
    #         (1, 4, numpy.array([0.3, 0.5, -0.9])),
    #         (2, 3, numpy.array([0.1, -0.3, -0.4])),
    #         (2, 4, numpy.array([-0.1, -0.3, -1.4])),
    #         (3, 4, numpy.array([0.1, -0.6, -2.4])),
    #         ]

    # table = [ [[1, 2, 4], [2, 3, 4]],
    #         [[3, 4, 4]] ]

    for tt in table:
        # example: tt = [[1,2,4], [2,3,4]]
        iatm_beg, jatm_beg, jatm_end = tt[0]
        iatm_end = tt[-1][0]
        table_tbfs = numpy.zeros(
                (iatm_end-iatm_beg+1, natom-iatm_beg+1, 3), numpy.float)

        # print(iatm_beg, iatm_end, jatm_beg, jatm_end)
        # print(table_tbfs.shape)

        for iatm, jatm, tbf in pair_to_tbfs:
            if not (iatm_beg <= iatm <= iatm_end): continue
            if not (jatm_beg <= jatm <= jatm_end): continue
            # # table_tbfs[iatm-iatm_beg, jatm-iatm_beg, :] = tbf[:]

            # print(iatm, jatm, iatm-iatm_beg, jatm-iatm_beg)
            table_tbfs[iatm-iatm_beg, jatm-iatm_beg] = tbf
            # print(table_tbs)
            # table_tbs[:] = tbf[:]

        yield table_tbfs

def gen_table_tbfs(pair_to_tbfs, table, natom):

    # pair_to_tbfs = [
    #         (1, 2, numpy.array([0.2, 0.3, 0.4])),
    #         (1, 3, numpy.array([0.5, 0.6, 0.9])),
    #         (1, 4, numpy.array([0.3, 0.5, -0.9])),
    #         (2, 3, numpy.array([0.1, -0.3, -0.4])),
    #         (2, 4, numpy.array([-0.1, -0.3, -1.4])),
    #         (3, 4, numpy.array([0.1, -0.6, -2.4])),
    #         ]

    # table = [ [[1, 2, 4], [2, 3, 4]],
    #         [[3, 4, 4]] ]

    for tt in table:
        # example: tt = [[1,2,4], [2,3,4]]
        iatm_beg, jatm_beg, jatm_end = tt[0]
        iatm_end = tt[-1][0]
        table_tbfs = numpy.zeros(
                (iatm_end-iatm_beg+1, natom-iatm_beg+1, 3), numpy.float)

        for t in tt:
            iatm__, jatm_beg, jatm_end = t

            for iatm, jatm, tbf in pair_to_tbfs:
                if iatm__ != iatm: continue
                if not (jatm_beg <= jatm <= jatm_end): continue

                table_tbfs[iatm-iatm_beg, jatm-iatm_beg] = tbf

        yield table_tbfs

def gen_pair_to_tbfs_with_table(pair_to_tbfs, table):

    # pair_to_tbfs = [
    #         (1, 2, numpy.array([0.2, 0.3, 0.4])),
    #         (1, 3, numpy.array([0.5, 0.6, 0.9])),
    #         (1, 4, numpy.array([0.3, 0.5, -0.9])),
    #         (2, 3, numpy.array([0.1, -0.3, -0.4])),
    #         (2, 4, numpy.array([-0.1, -0.3, -1.4])),
    #         (3, 4, numpy.array([0.1, -0.6, -2.4])),
    #         ]

    # table = [ [[1, 2, 4], [2, 3, 4]],
    #         [[3, 4, 4]] ]

    for tt in table:
        # example: tt = [[1,2,4], [2,3,4]]
        # iatm_beg, jatm_beg, jatm_end = tt[0]
        # iatm_end = tt[-1][0]

        for t in tt:
            iatm_beg, jatm_beg, jatm_end = t

            for iatm, jatm, tbf in pair_to_tbfs:
                if iatm_beg != iatm: continue
                if not (jatm_beg <= jatm <= jatm_end): continue

                yield (iatm, jatm, tbf)

def cal_pot_current(crd, iatm_to_itars, pair_to_tbfs, volumes):
    ntarget = max(iatm_to_itars)
    stress_tensor  = numpy.zeros([ntarget, 3, 3])
    for iatm, jatm, f_ij in pair_to_tbfs:
        itar = iatm_to_itars[iatm-1]
        jtar = iatm_to_itars[jatm-1]

        # for iatm
        r_i, r_j = crd[iatm-1], crd[jatm-1]
        if itar != 0:
            st_i = stress_tensor[itar-1]
            st_i += 0.5 * tensor(f_ij, r_i-r_j) / volumes[itar-1]

        # for jatm
        if jtar != 0:
            st_j = stress_tensor[jtar-1]
            st_j += 0.5 * tensor(-f_ij, r_j-r_i) / volumes[jtar-1]

    return stress_tensor


################################################################################
from current.current import StressTensor
class TestStressTensorWithTarget:

    def test_kinetic_of_3_particle(self):

        # preparing
        vel = numpy.array([
            [0.11, 0.08, -0.02],
            [0.33, 0.23, -0.11],
            [-0.13, 0.03, 0.25], ] )
        masses  = [1.08, 16.01, 1.08]
        volumes = [8.0, 20.0, 8.0]

        natom = len(vel)
        target_atoms = [1,3]
        iatm_to_itars = get_iatm_to_itars(target_atoms, natom)
        iatm_to_igrps = numpy.zeros([natom], numpy.int)

        target_vols = [ volumes[iatm-1] for iatm in target_atoms ]

        # answer
        answer = list(gen_kinetic_current(vel, masses,
                        target_vols, iatm_to_itars))

        # result
        scal = StressTensor(target_atoms, iatm_to_igrps)
        cur_atm, cur_inn, cur_out = scal.cal_kinetic(vel, masses, target_vols)
        result = cur_atm

        print(answer)
        print(result)
        npt.assert_almost_equal(result, answer)

    def test_bonded_of_2_particles(self):

        crd = numpy.array(
            [[1.0, 2.0, 3.0],
             [2.0, 1.0, 1.5]] )
        pair_to_tbfs = [(1, 2, numpy.array([0.2, 0.3, 0.4]))]
        volumes = [1.2, 1.4]
        natom = len(crd)
        target_atoms = [2]
        iatm_to_itars = get_iatm_to_itars(target_atoms, natom)
        iatm_to_igrps = [1,1]

        target_vols = [ volumes[iatm-1] for iatm in target_atoms ]

        # answer
        answer = cal_pot_current(crd, iatm_to_itars, pair_to_tbfs, target_vols)

        scal = StressTensor(target_atoms, iatm_to_igrps)
        tbfs, zero_tbfs = get_extended_tbfs(pair_to_tbfs, natom=natom, n_ext=10)
        values = scal.cal_bonded(crd, tbfs, target_vols, zero_tbfs)
        cur_atm, cur_inn, cur_out = values
        result = cur_atm

        # print(result)
        # print(answer)
        npt.assert_almost_equal(result, answer)

    def test_nonbonded_of_2_particles(self):
        crd = numpy.array(
            [[1.0, 2.0, 3.0],
             [2.0, 1.0, 1.5]] )
        pair_to_tbfs = [(1, 2, numpy.array([0.2, 0.3, 0.4]))]
        volumes = [1.2, 1.4]

        natom = len(crd)
        target_atoms = [2]
        iatm_to_itars = get_iatm_to_itars(target_atoms, natom)

        # group
        iatm_to_igrps = [1]*natom
        gvolumes = [1.0]

        # table
        table = [[[1, 2, 2]]]

        target_vols = [ volumes[iatm-1] for iatm in target_atoms ]

        # right answer
        new_pair_to_tbfs = gen_pair_to_tbfs_with_table(pair_to_tbfs, table)
        answer = cal_pot_current(crd, iatm_to_itars, pair_to_tbfs, target_vols)

        # result
        scal = StressTensor(target_atoms, iatm_to_igrps)
        table_tbfs = list(gen_table_tbfs(pair_to_tbfs, table, natom=natom))
        cur_atm, cur_inn, cur_out = scal.cal_nonbonded(
                                crd, table_tbfs, target_vols, table, gvolumes)
        result = cur_atm

        # assert
        npt.assert_almost_equal(result, answer)

        # print(result)
        # print(answer)

        # ok_(False)


    def test_nonbonded_of_4_particles(self):
        crd = numpy.array(
            [[1.0, 2.0,   3.0],
             [2.0, 1.0,   1.5],
             [2.5, -1.0, -1.0],
             [-0.2, -0.3, -0.6]] )
        pair_to_tbfs = [
                (1, 2, numpy.array([0.2, 0.3, 0.4])),
                (1, 3, numpy.array([0.5, 0.6, 0.9])),
                (1, 4, numpy.array([0.3, 0.5, -0.9])),
                (2, 3, numpy.array([0.1, -0.3, -0.4])),
                (2, 4, numpy.array([-0.1, -0.3, -1.4])),
                (3, 4, numpy.array([0.1, -0.6, -2.4])),
                ]
        volumes = [1.2, 1.4, 1.3, 1.6]

        natom = len(crd)
        target_atoms = [2,3,4]
        iatm_to_itars = get_iatm_to_itars(target_atoms, natom)

        # table
        table = [
            [[1, 2, 4], [2, 3, 3]],
            [[2, 4, 4], [3, 4, 4]] ]

        # group
        iatm_to_igrps = [1]*natom
        gvolumes = [1.0]

        target_vols = [ volumes[iatm-1] for iatm in target_atoms ]

        # right answer
        new_pair_to_tbfs = gen_pair_to_tbfs_with_table(pair_to_tbfs, table)
        answer = cal_pot_current(crd, iatm_to_itars, pair_to_tbfs, target_vols)

        # result
        scal = StressTensor(target_atoms, iatm_to_igrps)
        table_tbfs = list(gen_table_tbfs(pair_to_tbfs, table, natom=natom))
        cur_atm, cur_inn, cur_out = scal.cal_nonbonded(
                                crd, table_tbfs, target_vols, table, gvolumes)
        result = cur_atm

        # assert
        npt.assert_almost_equal(result, answer)

        # print(result)
        # print(answer)

        # ok_(False)



# class TestStressTensor:

#     # def test_1(self):
#     #     assert 1 == 1

#     def test_kinetic_of_1_particle(self):
#         scal = StressTensor()
#         vel = numpy.array([[0.11,0.08,-0.02]])
#         masses  = [1.08]
#         volumes = [1.2]

#         right_cur = gen_kinetic_current(vel, volumes, masses)
#         res  = scal.cal_kinetic(vel, volumes, masses)

#         for r, rc in zip(res , right_cur):
#             npt.assert_almost_equal(r, rc)

    # def test_kinetic_of_3_particles(self):
    #     scal = StressTensor()
    #     vel = numpy.array([[0.11,0.08,-0.02],
    #                        [0.21,-0.10,0.05],
    #                        [-0.05, 0.22, 0.15]])
    #     masses  = [1.08, 18.02, 1.08]
    #     volumes = [1.2, 1.4, 1.2]

    #     right_cur = gen_kinetic_current(vel, volumes, masses)
    #     res  = scal.cal_kinetic(vel, volumes, masses)

    #     for r, rc in zip(res , right_cur):
    #         npt.assert_almost_equal(r, rc)

    # def test_bonded_of_2_particles(self):

    #     crd = numpy.array(
    #         [[1.0, 2.0, 3.0],
    #          [2.0, 1.0, 1.5]] )
    #     pair_to_tbfs = [(1, 2, numpy.array([0.2, 0.3, 0.4]))]
    #     volumes = [1.2, 1.4]
    #     natom = len(crd)
    #     target_atoms = [1]
    #     iatm_to_itars = get_iatm_to_itars(target_atoms, natom)

    #     # answer
    #     answer = cal_pot_current(crd, iatm_to_itars, pair_to_tbfs, volumes)

    #     scal = StressTensor(target_atoms, iatm_to_igrps=[])
    #     tbfs, zero_tbfs = get_extended_tbfs(pair_to_tbfs, natom=natom, n_ext=10)
    #     values = scal.cal_bonded(crd, tbfs, volumes, zero_tbfs)
    #     cur_atm, cur_inn, cur_out = values
    #     result = cur_atm

    #     print(result)
    #     print(answer)
    #     npt.assert_almost_equal(result, answer)

    # def test_bonded_of_3_particles(self):
    #     crd = numpy.array(
    #         [[1.0, 2.0,   3.0],
    #          [2.0, 1.0,   1.5],
    #          [2.5, -1.0, -1.0]] )
    #     pair_to_tbfs = [
    #             (1, 2, numpy.array([0.2, 0.3, 0.4])),
    #             (1, 3, numpy.array([0.5, 0.6, 0.9])),
    #             (2, 3, numpy.array([0.1, -0.3, -0.4]))
    #             ]
    #     volumes = [1.2, 1.4, 1.2]

    #     right = cal_pot_current(crd, pair_to_tbfs, volumes)

    #     scal = StressTensor()
    #     tbfs, zero_tbfs = get_extended_tbfs(pair_to_tbfs, natom=3, n_ext=10)
    #     res = scal.cal_bonded(crd, tbfs, volumes, zero_tbfs)

    #     # npt.assert_almost_equal(res, right)
    #     npt.assert_almost_equal(res, right)

    # def test_bonded_of_4_particles(self):
    #     crd = numpy.array(
    #         [[1.0, 2.0,   3.0],
    #          [2.0, 1.0,   1.5],
    #          [2.5, -1.0, -1.0],
    #          [-0.2, -0.3, -0.6]] )
    #     pair_to_tbfs = [
    #             (1, 2, numpy.array([0.2, 0.3, 0.4])),
    #             (1, 3, numpy.array([0.5, 0.6, 0.9])),
    #             (1, 4, numpy.array([0.3, 0.5, -0.9])),
    #             (2, 3, numpy.array([0.1, -0.3, -0.4])),
    #             (2, 4, numpy.array([-0.1, -0.3, -1.4])),
    #             (3, 4, numpy.array([0.1, -0.6, -2.4])),
    #             ]
    #     volumes = [1.2, 1.4, 1.2, 1.6]

    #     right = cal_pot_current(crd, pair_to_tbfs, volumes)

    #     scal = StressTensor()
    #     tbfs, zero_tbfs = get_extended_tbfs(pair_to_tbfs, natom=4, n_ext=10)
    #     res = scal.cal_bonded(crd, tbfs, volumes, zero_tbfs)

    #     npt.assert_almost_equal(res, right)

    # def test_nonbonded_of_2_particles(self):
    #     crd = numpy.array(
    #         [[1.0, 2.0, 3.0],
    #          [2.0, 1.0, 1.5]] )
    #     pair_to_tbfs = [(1, 2, numpy.array([0.2, 0.3, 0.4]))]
    #     volumes = [1.2, 1.4]

    #     table = [[[1, 2, 2]]]

    #     # right answer
    #     new_pair_to_tbfs = gen_pair_to_tbfs_with_table(pair_to_tbfs, table)
    #     right = cal_pot_current(crd, new_pair_to_tbfs, volumes)

    #     # result
    #     scal = StressTensor()
    #     table_tbfs = list(gen_table_tbfs(pair_to_tbfs, table, natom=2))
    #     res = scal.cal_nonbonded(crd, table_tbfs, volumes, table)

    #     # assert
    #     npt.assert_almost_equal(res, right)

    # def test_nonbonded_of_4_particles(self):
    #     crd = numpy.array(
    #         [[1.0, 2.0,   3.0],
    #          [2.0, 1.0,   1.5],
    #          [2.5, -1.0, -1.0],
    #          [-0.2, -0.3, -0.6]] )
    #     pair_to_tbfs = [
    #             (1, 2, numpy.array([0.2, 0.3, 0.4])),
    #             (1, 3, numpy.array([0.5, 0.6, 0.9])),
    #             (1, 4, numpy.array([0.3, 0.5, -0.9])),
    #             (2, 3, numpy.array([0.1, -0.3, -0.4])),
    #             (2, 4, numpy.array([-0.1, -0.3, -1.4])),
    #             (3, 4, numpy.array([0.1, -0.6, -2.4])),
    #             ]
    #     volumes = [1.2, 1.4, 1.2, 1.6]

    #     table = [
    #         [[1, 2, 4], [2, 3, 3]],
    #         [[2, 4, 4], [3, 4, 4]] ]

    #     # right answer
    #     new_pair_to_tbfs = list(
    #             gen_pair_to_tbfs_with_table(pair_to_tbfs, table))
    #     right = cal_pot_current(crd, new_pair_to_tbfs, volumes)
    #     print(numpy.array(right))
    #     print(new_pair_to_tbfs)
    #     print('*'*80)

    #     # result
    #     scal = StressTensor()
    #     table_tbfs = list(gen_table_tbfs(pair_to_tbfs, table, natom=4))
    #     res = scal.cal_nonbonded(crd, table_tbfs, volumes, table)
    #     print(res)

    #     # assert
    #     npt.assert_almost_equal(res, right)

    # def test_nonbonded_of_4_particles_table2(self):
    #     crd = numpy.array(
    #         [[1.0, 2.0,   3.0],
    #          [2.0, 1.0,   1.5],
    #          [2.5, -1.0, -1.0],
    #          [-0.2, -0.3, -0.6]] )
    #     pair_to_tbfs = [
    #             (1, 2, numpy.array([0.2, 0.3, 0.4])),
    #             (1, 3, numpy.array([0.5, 0.6, 0.9])),
    #             (1, 4, numpy.array([0.3, 0.5, -0.9])),
    #             (2, 3, numpy.array([0.1, -0.3, -0.4])),
    #             (2, 4, numpy.array([-0.1, -0.3, -1.4])),
    #             (3, 4, numpy.array([0.1, -0.6, -2.4])),
    #             ]
    #     volumes = [1.2, 1.4, 1.2, 1.6]

    #     table = [
    #         [[1, 2, 4]],
    #         [[2, 3, 3]],
    #         [[2, 4, 4], [3, 4, 4]] ]

    #     # right answer
    #     new_pair_to_tbfs = gen_pair_to_tbfs_with_table(pair_to_tbfs, table)
    #     right = cal_pot_current(crd, new_pair_to_tbfs, volumes)

    #     # result
    #     scal = StressTensor()
    #     table_tbfs = list(gen_table_tbfs(pair_to_tbfs, table, natom=4))
    #     res = scal.cal_nonbonded(crd, table_tbfs, volumes, table)

    #     # assert
    #     npt.assert_almost_equal(res, right)

    # def test_nonbonded_of_4_particles_table3(self):
    #     crd = numpy.array(
    #         [[1.0, 2.0,   3.0],
    #          [2.0, 1.0,   1.5],
    #          [2.5, -1.0, -1.0],
    #          [-0.2, -0.3, -0.6]] )
    #     pair_to_tbfs = [
    #             (1, 2, numpy.array([0.2, 0.3, 0.4])),
    #             (1, 3, numpy.array([0.5, 0.6, 0.9])),
    #             (1, 4, numpy.array([0.3, 0.5, -0.9])),
    #             (2, 3, numpy.array([0.1, -0.3, -0.4])),
    #             (2, 4, numpy.array([-0.1, -0.3, -1.4])),
    #             (3, 4, numpy.array([0.1, -0.6, -2.4])),
    #             ]
    #     volumes = [1.2, 1.4, 1.2, 1.6]

    #     table = [
    #         [[1, 4, 4]],
    #         [[2, 3, 3], [3, 4, 4]]]

    #     # right answer
    #     new_pair_to_tbfs = gen_pair_to_tbfs_with_table(pair_to_tbfs, table)
    #     right = cal_pot_current(crd, new_pair_to_tbfs, volumes)

    #     # result
    #     scal = StressTensor()
    #     table_tbfs = list(gen_table_tbfs(pair_to_tbfs, table, natom=4))
    #     res = scal.cal_nonbonded(crd, table_tbfs, volumes, table)

    #     # assert
    #     npt.assert_almost_equal(res, right)

