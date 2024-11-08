from nose.tools import *
import numpy.testing as npt

import os, sys
import numpy as np
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)

# module to test
from volume import calvolume as volmod

################################################################################
class TestVoronoiVolume:

    def setup(self):

        amber_dir = os.path.abspath( os.path.join(
            os.path.dirname(__file__), '..', '..', 'parser', 'amber'))
        if amber_dir not in sys.path:
            sys.path.insert(0, amber_dir)
        from topology import TopologyParser, Format2Amber99Converter
        from trajectory import CoordinateParser

        traj_fn = os.path.join(amber_dir, 'test/ala3-woct.mdcrd.gz')
        tpl_fn  = os.path.join(amber_dir, 'test/ala3-woct.prmtop')

        # prepare topology
        tpl = TopologyParser(tpl_fn)
        tpl.parse()
        conv = Format2Amber99Converter(tpl)
        conv.convert()

        # prepare trajectory
        traj_parser = CoordinateParser(traj_fn, 5844)
        crd, box = next(traj_parser)

        # target_atoms = [1,5,7,11,12,13,15,17,21,22,23,25,27,31,32,33]
        target_atoms = list(range(1, 33+1))
        grp_to_atoms = [
                ('ALA_001', list(range(1, 12+1))),
                ('ALA_002', list(range(13, 22+1))),
                ('ALA_003', list(range(23, 33+1))), ]
        vsetting = volmod.VolumeSetting(topology=conv, traj_parser=traj_parser,
                target_atoms=target_atoms, grp_to_atoms=grp_to_atoms )

        return vsetting, crd

    def setup_with_united(self):

        amber_dir = os.path.abspath( os.path.join(
            os.path.dirname(__file__), '..', '..', 'parser', 'amber'))
        if amber_dir not in sys.path:
            sys.path.insert(0, amber_dir)
        from topology import TopologyParser, Format2Amber99Converter
        from trajectory import CoordinateParser

        traj_fn = os.path.join(amber_dir, 'test/ala3-woct.mdcrd.gz')
        tpl_fn  = os.path.join(amber_dir, 'test/ala3-woct.prmtop')

        # prepare topology
        tpl = TopologyParser(tpl_fn)
        tpl.parse()
        conv = Format2Amber99Converter(tpl)
        conv.convert()

        # prepare trajectory
        traj_parser = CoordinateParser(traj_fn, 5844)
        crd, box = next(traj_parser)

        # target_atoms = [1,5,7,11,12,13,15,17,21,22,23,25,27,31,32,33]
        target_atoms = list(range(1, 33+1))
        grp_to_atoms = [
                ('00001_N', [1, 2, 3, 4]),
                ('00005_CA', [5, 6]),
                ('00007_CB', [7, 8, 9, 10]),
                ('00011_C', [11]),
                ('00012_O', [12]),
                ('00013_N', [13, 14]),
                ('00015_CA', [15, 16]),
                ('00017_CB', [17, 18, 19, 20]),
                ('00021_C', [21]),
                ('00022_O', [22]),
                ('00023_N', [23, 24]),
                ('00025_CA', [25, 26]),
                ('00027_CB', [27, 28, 29, 30]),
                ('00031_C', [31]),
                ('00032_O', [32]),
                ('00033_OXT', [33]) ]

        vsetting = volmod.VolumeSetting(topology=conv, traj_parser=traj_parser,
                target_atoms=target_atoms, grp_to_atoms=grp_to_atoms )

        return vsetting, crd

    def test_fast_and_debug(self): 
        # preparation
        vsetting, crd = self.setup()
        vsetting.voronoi_no_hydrogen = False
        vsetting.voronoi_cutoff = 6

        vcal = volmod.VoronoiVolumeCalculator(vsetting)
        
        # check the fast method and the vervose method
        radii, volumes     = vcal.cal_voronoi(crd)
        radii_d, volumes_d = vcal.cal_voronoi_with_debug(crd)

        # npt.assert_almost_equal(volumes, volumes_d)
        npt.assert_almost_equal(volumes, volumes_d)

    def test_hydrogen(self):

        # preparation
        vsetting, crd = self.setup()
        vsetting.voronoi_no_hydrogen = False
        vsetting.voronoi_cutoff = 6

        vcal = volmod.VoronoiVolumeCalculator(vsetting)

        # answer value
        answer_atom = np.array(
                [2.53095618, 10.57237302, 8.77409149, 10.28907996, 4.39428955,
                 9.56204114, 2.86505444, 14.06828119, 13.29667321, 14.26201939,
                 6.91872261, 9.26958893, 5.98605445, 9.55967883, 4.38770271,
                 16.75400934, 2.77741587, 12.69864899, 18.34138455, 11.70118055,
                 7.42848285, 10.79324531, 6.64302911, 12.19472874, 4.76973015,
                 13.0583563, 3.00119366, 12.10096405,  10.79443561, 14.06985516,
                 7.16893742, 6.12357164, 7.19708035] )

        answer_group = np.array([106.80317111, 100.42780345, 97.12188219])
        
        # check volumes
        volumes  = vcal.get_volume(crd)
        gvolumes = vcal.get_gvolume(crd)

        npt.assert_almost_equal(volumes,  answer_atom)
        npt.assert_almost_equal(gvolumes, answer_group)

    def test_no_hydrogen(self):

        # preparation
        vsetting, crd = self.setup()
        vsetting.voronoi_no_hydrogen = True
        vsetting.voronoi_cutoff = 6

        vcal = volmod.VoronoiVolumeCalculator(vsetting)

        # answer value
        answer_atom = np.array(
            [20.22022312 , 8.          , 8.          , 8. , 10.98203136 , 
            8.           , 30.56816177 , 8.          , 8. , 8.          , 
            8.12470536   , 18.52584001 , 12.99162867 , 8. , 15.54317338 , 
            8.           , 33.70637798 , 8.          , 8. , 8.          , 
            8.39924737   , 17.52766371 , 15.78170899 , 8. , 13.83811941 , 
            8.           , 27.97682307 , 8.          , 8. , 8.          , 
            8.35937656   , 16.07189061 , 16.88228382] )

        answer_group = np.array([144.42096161, 128.16809111, 138.91020245])
        
        # check volumes
        volumes  = vcal.get_volume(crd)
        gvolumes = vcal.get_gvolume(crd)

        npt.assert_almost_equal(volumes,  answer_atom)
        npt.assert_almost_equal(gvolumes, answer_group)

    def test_no_hydrogen_with_united(self):

        # preparation
        vsetting, crd = self.setup_with_united()
        vsetting.voronoi_no_hydrogen = True
        vsetting.group_method = 'united'
        vsetting.voronoi_cutoff = 6

        vcal = volmod.VoronoiVolumeCalculator(vsetting)

        # answer value
        answer_atom = np.array(
            [20.22022312 , 8.          , 8.          , 8. , 10.98203136 , 
            8.           , 30.56816177 , 8.          , 8. , 8.          , 
            8.12470536   , 18.52584001 , 12.99162867 , 8. , 15.54317338 , 
            8.           , 33.70637798 , 8.          , 8. , 8.          , 
            8.39924737   , 17.52766371 , 15.78170899 , 8. , 13.83811941 , 
            8.           , 27.97682307 , 8.          , 8. , 8.          , 
            8.35937656   , 16.07189061 , 16.88228382] )

        answer_group = np.array(
            [20.22022312 , 10.98203136 , 30.56816177 ,  8.12470536 ,
             18.52584001 , 12.99162867 , 15.54317338 , 33.70637798 ,
              8.39924737 , 17.52766371 , 15.78170899 , 13.83811941 , 
             27.97682307 , 8.35937656  , 16.07189061 , 16.88228382] )
        
        # check volumes
        volumes  = vcal.get_volume(crd)
        gvolumes = vcal.get_gvolume(crd)

        npt.assert_almost_equal(volumes,  answer_atom)
        npt.assert_almost_equal(gvolumes, answer_group)


class TestVoronoiVolumeWithSolvation:

    def setup(self):

        amber_dir = os.path.abspath( os.path.join(
            os.path.dirname(__file__), '..', '..', 'parser', 'amber'))
        if amber_dir not in sys.path:
            sys.path.insert(0, amber_dir)
        from topology import TopologyParser, Format2Amber99Converter
        from trajectory import CoordinateParser

        traj_fn = os.path.join(amber_dir, 'test/ala3.mdcrd.gz')
        tpl_fn  = os.path.join(amber_dir, 'test/ala3.prmtop')

        # prepare topology
        tpl = TopologyParser(tpl_fn)
        tpl.parse()
        conv = Format2Amber99Converter(tpl)
        conv.convert()

        # prepare trajectory
        traj_parser = CoordinateParser(traj_fn, 33)
        crd, box = next(traj_parser)

        # target_atoms = [1,5,7,11,12,13,15,17,21,22,23,25,27,31,32,33]
        target_atoms = list(range(1, 33+1))
        grp_to_atoms = [
                ('ALA_001', list(range(1, 12+1))),
                ('ALA_002', list(range(13, 22+1))),
                ('ALA_003', list(range(23, 33+1))), ]
        vsetting = volmod.VolumeSetting(topology=conv, traj_parser=traj_parser,
                target_atoms=target_atoms, grp_to_atoms=grp_to_atoms )

        vsetting.voronoi_solvation = 'RANDOM20'
        vsetting.voronoi_cutoff = 6

        return vsetting, crd

    def test_fast_and_debug(self): 
        # preparation
        vsetting, crd = self.setup()
        vsetting.voronoi_no_hydrogen = False

        vcal = volmod.VoronoiVolumeCalculatorWithSolvation(vsetting)
        
        # check the fast method and the vervose method
        radii, volumes     = vcal.cal_voronoi(crd)
        radii_d, volumes_d = vcal.cal_voronoi_with_debug(crd)

        # npt.assert_almost_equal(volumes, volumes_d)
        npt.assert_almost_equal(volumes, volumes_d)

    def test_hydrogen(self):

        # preparation
        vsetting, crd = self.setup()
        vsetting.voronoi_no_hydrogen = False

        vcal = volmod.VoronoiVolumeCalculatorWithSolvation(vsetting)

        # answer value
        answer_atom = np.array(
                [2.47223895, 19.14232321, 16.21803614,  7.92137772,  4.22816654,
                19.82479296,  3.13903811, 12.11502524, 14.1649559,  15.90249363,
                 5.93714983, 13.81966351,  4.72193302,  5.39133486,  3.58498268,
                 5.56483001,  2.51259541,  4.1971505 ,  5.93018255,  4.54570011,
                 5.86307023, 13.2699581 ,  5.28416826,  5.84329651,  4.512589,
                15.45197065,  2.80961053, 14.85935785, 18.02848079,  8.58282678,
                 7.61211716, 12.44204141, 18.46705607] )

        answer_group = np.array([134.885261726, 55.581737472, 113.893515011])

        # check volumes
        volumes  = vcal.get_volume(crd)
        gvolumes = vcal.get_gvolume(crd)

        npt.assert_almost_equal(volumes,  answer_atom)
        npt.assert_almost_equal(gvolumes, answer_group)

    def test_no_hydrogen(self):

        # preparation
        vsetting, crd = self.setup()
        vsetting.voronoi_no_hydrogen = True

        vcal = volmod.VoronoiVolumeCalculatorWithSolvation(vsetting)

        # answer value
        answer_atom = np.array(
            [27.4864788  , 8.          , 8.          , 8. , 14.37028708 ,
              8.         , 26.71430444 , 8.          , 8. ,  8.         ,
              7.9492852  , 20.02736218 , 9.12244069  , 8. ,  6.2135348  ,
              8.         , 11.20737756 , 8.          , 8. ,  8.         ,
              7.1022755  , 20.66406397 , 9.74530577  , 8. , 12.24821661 ,
              8.         , 33.86800002 , 8.          , 8. ,  8.         ,
              8.48899055 , 20.38073956 , 19.72677428] )

        answer_group = np.array([152.54771770, 94.309692516, 144.458026778])
        # check volumes
        volumes  = vcal.get_volume(crd)
        gvolumes = vcal.get_gvolume(crd)

        npt.assert_almost_equal(volumes,  answer_atom)
        npt.assert_almost_equal(gvolumes, answer_group)

    def test_solvation(self):
        """Ok if the system coordinate is solvated correctly."""

        # preparation
        vsetting, crd = self.setup()

        vcal = volmod.VoronoiVolumeCalculatorWithSolvation(vsetting)

        natom = vsetting.natom
        nres = max(vsetting.rids)

        # perform
        new_crd = vcal.solvate_crd(crd)
        natom_shell = len(new_crd) - natom
        nres_shell = natom_shell/3

        len2 = 1.10**2

        def dist2(pos1, pos2):
            r_12 = pos1 - pos2
            return r_12[0]*r_12[0] + r_12[1]*r_12[1] + r_12[2]*r_12[2]

        is_water_positions = []
        # for ires in range(nres+1, nres+nres_shell+1):
        #     iatm_o_1 = natom + 3*(ires-1)
        for ires in range(nres_shell):
            iatm_o_1 = natom + 3*ires
            
            pos_o  = new_crd[iatm_o_1]
            pos_1 = new_crd[iatm_o_1+1]
            pos_2 = new_crd[iatm_o_1+2]

            len2_o_1 = dist2(pos_o, pos_1)
            len2_o_2 = dist2(pos_o, pos_2)
            len2_1_2 = dist2(pos_1, pos_2)

            if (len2_o_1 < len2) and (len2_o_2 < len2):
                is_water_positions.append( True )
            else:
                is_water_positions.append( False )

        # answer value
        answer_value = [ True for ires in range(nres_shell) ]
        
        # check
        eq_(natom_shell%3, 0)
        eq_(is_water_positions, answer_value)


################################################################################
class TestVolume1:

    def setup(self):

        # target_atoms = [1,5,7,11,12,13,15,17,21,22,23,25,27,31,32,33]
        natom = 33
        target_atoms = list(range(1, natom+1))
        grp_to_atoms = [
                ('ALA_001', list(range(1, 12+1))),
                ('ALA_002', list(range(13, 22+1))),
                ('ALA_003', list(range(23, natom+1))), ]

        vsetting = volmod.VolumeSetting(
                target_atoms=target_atoms, grp_to_atoms=grp_to_atoms )

        vsetting.natom = natom

        return vsetting

    def test_volume1(self): 

        # preparation
        vsetting = self.setup()
        vcal = volmod.Volume1(vsetting)
        crd = []

        # answer value
        answer_atom = np.ones([vsetting.natom])
        answer_group = np.ones([len(vsetting.grp_to_atoms)])
        
        # perform
        volumes  = vcal.get_volume()
        gvolumes = vcal.get_gvolume()

        # check
        npt.assert_almost_equal(volumes,  answer_atom)
        npt.assert_almost_equal(gvolumes, answer_group)

################################################################################
class TestVDWVolume:

    def setup(self):

        # target_atoms = [1,5,7,11,12,13,15,17,21,22,23,25,27,31,32,33]
        natom = 33
        target_atoms = list(range(1, natom+1))
        grp_to_atoms = [
                ('ALA_001', list(range(1, 12+1))),
                ('ALA_002', list(range(13, 22+1))),
                ('ALA_003', list(range(23, natom+1))), ]

        vdw_radii = [ 1.86, 0.6, 0.6, 0.6, 1.9, 0.6, 1.9, 0.6, 0.6, 0.6,
                       1.9, 1.8, 1.86, 0.6, 1.9, 0.6, 1.9, 0.6, 0.6, 0.6,
                       1.9, 1.8, 1.86, 0.6, 1.9, 0.6, 1.9, 0.6, 0.6, 0.6,
                       1.9, 1.8, 1.8 ]

        vsetting = volmod.VolumeSetting(
                target_atoms=target_atoms, grp_to_atoms=grp_to_atoms )

        vsetting.natom = natom
        vsetting.vdw_radii = vdw_radii

        return vsetting

    def test_vdw(self): 

        # preparation
        vsetting = self.setup()
        vcal = volmod.VDWVolumeCalculator(vsetting)

        # answer value
        answer_atom = np.array(
            [ 26.95426178,  0.90477868,  0.90477868,  0.90477868, 28.73091201,
               0.90477868, 28.73091201,  0.90477868,  0.90477868,  0.90477868,
              28.73091201, 24.42902447, 26.95426178,  0.90477868, 28.73091201,
               0.90477868, 28.73091201,  0.90477868,  0.90477868,  0.90477868,
              28.73091201, 24.42902447, 26.95426178,  0.90477868, 28.73091201,
               0.90477868, 28.73091201,  0.90477868,  0.90477868,  0.90477868,
              28.73091201, 24.42902447, 24.42902447] )

        answer_group = np.array(
            [143.90947308985176, 142.09991572138404, 166.52894019569828 ] )
        
        # perform
        volumes  = vcal.get_volume()
        gvolumes = vcal.get_gvolume()

        # check
        npt.assert_almost_equal(volumes,  answer_atom)
        npt.assert_almost_equal(gvolumes, answer_group)

################################################################################
class TestOuterVolume:

    def setup(self):

        # preparation
        natom = 33
        target_atoms = list(range(1, natom+1))
        vsetting = volmod.VolumeSetting()
        vsetting.natom = natom
        vsetting.target_atoms = target_atoms

        return vsetting

    def test_atomic_trajectory_only(self):

        # preparation
        vsetting = self.setup()
        vsetting.atomic_trajectory_file = os.path.join(
                os.path.dirname(__file__), 'dummy_vol.traj')
        vsetting.grp_to_atoms = [ 
                ('ALA_001', list(range(1, 12+1))),
                ('ALA_002', list(range(13, 22+1))),
                ('ALA_003', list(range(23, 33+1))), ]

        vcal = volmod.OuterVolumeFetcher(vsetting)
        
        # ** 1th step **
        # answer value
        answer_atom = np.array(
            [ 18.5, 10.1, 10.1, 10.1, 21.3, 8.2, 20.2, 8.2, 8.2, 8.2, 
              19. , 22. , 15.7,  7.7, 21.6, 8.4, 19.1, 7. , 7. , 7. , 
              19.85, 23., 17.4, 9.9 , 21.5, 8.4, 18.1, 8.8, 8.8, 8.8 , 
              19.7, 21.22, 21.43
        ])

        answer_group = np.array( [164.1, 136.35, 164.05] )

        # perform
        volumes  = vcal.get_volume()
        gvolumes = vcal.get_gvolume()

        # check
        npt.assert_almost_equal(volumes,  answer_atom)
        npt.assert_almost_equal(gvolumes, answer_group)

        # ** 2th step **
        # answer value
        answer_atom = np.array(
            [19.50, 9.10, 9.10, 9.10, 22.30, 8.20, 21.20, 7.20, 7.20, 7.20,
             20.00, 21.00, 18.70, 8.70, 22.60, 9.40, 18.10, 5.00, 5.00, 5.00,
             20.85, 24.00, 18.40, 9.70, 20.50, 9.40, 15.10, 7.80, 7.80, 7.80,
             18.70, 20.22, 20.43] )

        answer_group = np.array( [ 161.1 ,  137.35,  155.85] )

        # perform
        volumes  = vcal.get_volume()
        gvolumes = vcal.get_gvolume()

        # check
        npt.assert_almost_equal(volumes,  answer_atom)
        npt.assert_almost_equal(gvolumes, answer_group)

        # ** 3rd step **
        # answer value
        answer_atom = None
        answer_group = None

        # perform
        volumes  = vcal.get_volume()
        gvolumes = vcal.get_gvolume()

        # check
        eq_(volumes,  answer_atom)
        eq_(gvolumes, answer_group)


    def test_atomic_and_group(self):

        # preparation
        vsetting = self.setup()
        vsetting.atomic_trajectory_file = os.path.join(
                os.path.dirname(__file__), 'dummy_vol.traj')
        vsetting.group_trajectory_file = os.path.join(
                os.path.dirname(__file__), 'dummy_gvol.traj')
        vsetting.grp_to_atoms = [ 
               ('N_001'   , [1,2,3,4]     ) ,
               ('CA_005'  , [5,6]         ) ,
               ('CB_007'  , [7,8,9,10]    ) ,
               ('C_011'   , [11]          ) ,
               ('O_012'   , [12]          ) ,
               ('N_013'   , [13,14]       ) ,
               ('CA_015'  , [15,16]       ) ,
               ('CB_017'  , [17,18,19,20] ) ,
               ('C_021'   , [21]          ) ,
               ('O_022'   , [22]          ) ,
               ('N_023'   , [23,24]       ) ,
               ('CA_025'  , [25,26]       ) ,
               ('CB_027'  , [27,28,29,30] ) ,
               ('C_031'   , [31]          ) ,
               ('O_032'   , [32]          ) ,
               ('OXT_033' , [33]          ) 
        ]

        vcal = volmod.OuterVolumeFetcher(vsetting)
        
        # ** 1th step **
        # answer value
        answer_atom = np.array(
            [ 18.5, 10.1, 10.1, 10.1, 21.3, 8.2, 20.2, 8.2, 8.2, 8.2, 
              19. , 22. , 15.7,  7.7, 21.6, 8.4, 19.1, 7. , 7. , 7. , 
              19.85, 23., 17.4, 9.9 , 21.5, 8.4, 18.1, 8.8, 8.8, 8.8 , 
              19.7, 21.22, 21.43
        ])

        answer_group = np.array( 
                [18.50, 21.30, 20.20, 19.00, 22.00, 15.70, 21.60, 19.10, 19.85,
                 23.00, 17.40, 21.50, 18.10, 19.70, 21.22, 21.43] )

        # perform
        volumes  = vcal.get_volume()
        gvolumes = vcal.get_gvolume()

        # check
        npt.assert_almost_equal(volumes,  answer_atom)
        npt.assert_almost_equal(gvolumes, answer_group)

        # ** 2th step **
        # answer value
        answer_atom = np.array(
            [19.50, 9.10, 9.10, 9.10, 22.30, 8.20, 21.20, 7.20, 7.20, 7.20,
             20.00, 21.00, 18.70, 8.70, 22.60, 9.40, 18.10, 5.00, 5.00, 5.00,
             20.85, 24.00, 18.40, 9.70, 20.50, 9.40, 15.10, 7.80, 7.80, 7.80,
             18.70, 20.22, 20.43] )

        answer_group = np.array( 
                [18.50, 20.30, 19.20, 18.00, 20.00, 18.70, 20.60, 18.10, 19.85,
                 22.00, 18.40, 20.50, 19.10, 18.70, 22.22, 23.43] )

        # perform
        volumes  = vcal.get_volume()
        gvolumes = vcal.get_gvolume()

        # check
        npt.assert_almost_equal(volumes,  answer_atom)
        npt.assert_almost_equal(gvolumes, answer_group)

        # ** 3rd step **
        # answer value
        answer_atom = None
        answer_group = None

        # perform
        volumes  = vcal.get_volume()
        gvolumes = vcal.get_gvolume()

        # check
        eq_(volumes,  answer_atom)
        eq_(gvolumes, answer_group)

class TestVolumeBase:

    def test_setting(self):

        # preparation
        parameters = dict(
            natom = 50,
            elems = ['H', 'N', 'C'],
            grp_to_atoms = [
                    ('ALA_001', [1,2,3,4,5]),
                    ('ALA_002', [6,7,8,9,10]), ],
            target_atoms = [2,4,8,9],

            vdw_radii = [1.08, 1.80, 1.96],

            smve_rmax = 3.0,
            smve_dr   = 0.05,
            smve_interval = 2,
            smve_increment = 5,

            atomic_trajectory_file = 'volumes.traj',
            group_trajectory_file = 'gvolumes.traj',

            voronoi_cutoff       = 8.0,
            voronoi_no_hydrogen  = True,
            voronoi_solvation    = True,
            voronoi_probe_length = 2.6,
        )

        # set parameters
        vsetting = volmod.VolumeSetting()
        for key, value in list(parameters.items()):
            setattr(vsetting, key, value)

        # check
        for key, value in list(parameters.items()):
            print(key, getattr(vsetting, key), value)
            eq_(getattr(vsetting, key), value)

class TestSMVEVolume:

    # def test_smve(self):
    #     vsetting = volmod.VolumeSetting()
    #     vsetting.natom = 33
    #     vsetting.vdw_radii = None

    #     topdir = os.path.abspath( os.path.join(
    #         os.path.dirname(__file__), '..', '..', 'parser', 'amber'))
    #     if topdir not in sys.path:
    #         sys.path.insert(0, topdir)
    #     import trajectory
        
    #     traj_fn = os.path.join(topdir, 'test/ala3.mdcrd.gz')
    #     traj_parser = trajectory.CoordinateParser(traj_fn, vsetting.natom)
    #     vsetting.traj_parser = traj_parser
    #     vcal = volmod.SMVEVolumeCalculator(vsetting)

    #     parser2 = trajectory.CoordinateParser(traj_fn, vsetting.natom)
    #     crd, box = parser2.next()
    #     # crd, box = parser2.next()
    #     volumes = vcal.get_volume(crd)

    # def test_smve():
    #     vsetting = VolumeSetting()
    #     vsetting.natom = 1799
    #     vsetting.vdw_radii = None

    #     topdir = os.path.abspath( os.path.join(
    #         os.path.dirname(__file__), '..', 'parser', 'amber'))
    #     if topdir not in sys.path:
    #         sys.path.insert(0, topdir)
    #     import trajectory
        
    #     traj_fn = '../parser/amber/test/sam-nwat.mdcrd'
    #     traj_parser = trajectory.CoordinateParser(traj_fn, vsetting.natom)
    #     vsetting.traj_parser = traj_parser
    #     vcal = SMVEVolumeCalculator(vsetting)

    #     parser2 = trajectory.CoordinateParser(traj_fn, vsetting.natom)
    #     crd, box = parser2.next()
    #     # crd, box = parser2.next()
    #     volumes = vcal.get_volume(crd)

    pass

