from nose.tools import *
import numpy.testing as npt

import os, sys
import numpy
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)

# module to test
import group

################################################################################
class TestVolume:

    def test_with_target_atoms(self):

        # preparation
        natom = 33
        target_atoms = range(7,26+1)
        iatm_to_itars = group.get_iatm_to_itars(target_atoms, natom)

        res_info = dict(
            ids = 12*[1] + 10*[2] + 11*[3],
            names = 12*['ALA'] + 10*['GLU'] + 11*['LYS'],
        )

        atom_info = dict( 
            names = ['N','H1','H2','H3','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','N','H','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','N','H','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','OXT'], )

        # answer value
        answer1 = numpy.array(
                [ 0,  0,  0,  0,  0,  0,  1,  2,  3,  4,
                  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
                  15, 16, 17, 18, 19, 20,  0,  0,  0,  0,
                  0,  0,  0])

        answer2 = [
                ('00001_ALA', [7, 8, 9, 10, 11, 12]),
                ('00002_GLU', [13, 14, 15, 16, 17, 18, 19, 20, 21, 22]),
                ('00003_LYS', [23, 24, 25, 26]) ]

        # test
        res_pairs = list( group.gen_group_pairs_with_target(
                group.gen_residue_group(res_info, atom_info), iatm_to_itars ) )

        # check
        npt.assert_equal( iatm_to_itars, answer1)
        eq_( res_pairs, answer2)

    def test_file_group(self):

        # preparation
        filepath = os.path.join(topdir, 'tests/group.cfg')
        parser = group.GroupParser(filepath)

        # answer value
        answer = [
                ('00001_ALA', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]),
                ('00002_GLU', [13, 14, 15, 16, 17, 18, 19, 20, 21, 22]),
                ('00003_VAL', [23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33]),
        ]

        # test
        group_pairs = parser.gen_group_atoms_pairs()

        # check
        eq_( list(group_pairs), answer )

    def test_residue_group(self):

        # preparation
        res_info = dict(
            ids = 12*[1] + 10*[2] + 11*[3],
            names = 12*['ALA'] + 10*['ALA'] + 11*['ALA'],
        )

        atom_info = dict( 
            names = ['N','H1','H2','H3','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','N','H','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','N','H','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','OXT'], )

        # answer value
        answer = [
                ('00001_ALA', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]),
                ('00002_ALA', [13, 14, 15, 16, 17, 18, 19, 20, 21, 22]),
                ('00003_ALA', [23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33]),
        ]

        # test
        res_pairs = group.gen_residue_group(res_info, atom_info)

        # check
        eq_( list(res_pairs), answer )

    def test_united_group_1(self):

        # preparation
        atom_info = dict(
            names = ['N','H1','H2','H3','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','N','H','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','N','H','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','OXT'],

            ids   = range(1,33+1),

            elems = ['N','H','H','H','C','H','C','H','H','H',
                     'C','O','N','H','C','H','C','H','H','H',
                     'C','O','N','H','C','H','C','H','H','H',
                     'C','O','O'],
        )

        # answer value
        answer = [ ('00001_N', [1, 2, 3, 4]),
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
                   ('00033_OXT', [33])
        ]

        # test
        group_atoms_pairs = group.gen_united_group(atom_info)

        # check
        eq_( list(group_atoms_pairs), answer )


    def test_united_group_2(self):

        # preparation
        atom_info = dict(
            names = ['H1','N','H2','H3','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','N','H','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','N','H','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','OXT'],

            ids   = range(1,33+1),

            elems = ['H','N','H','H','C','H','C','H','H','H',
                     'C','O','N','H','C','H','C','H','H','H',
                     'C','O','N','H','C','H','C','H','H','H',
                     'C','O','O'],
        )

        # answer value
        answer = [ ('00002_N', [1, 2, 3, 4]),
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
                   ('00033_OXT', [33])
        ]

        # test
        group_atoms_pairs = group.gen_united_group(atom_info)

        # check
        eq_( list(group_atoms_pairs), answer )

    def test_united_group_3(self):

        # preparation
        atom_info = dict(
            names = ['H1','H2','N','H3','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','N','H','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','N','H','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','OXT'],

            ids   = range(1,33+1),

            elems = ['H','H','N','H','C','H','C','H','H','H',
                     'C','O','N','H','C','H','C','H','H','H',
                     'C','O','N','H','C','H','C','H','H','H',
                     'C','O','O'],
        )

        # answer value
        answer = [ ('00003_N', [1, 2, 3, 4]),
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
                   ('00033_OXT', [33])
        ]

        # test
        group_atoms_pairs = group.gen_united_group(atom_info)

        # check
        eq_( list(group_atoms_pairs), answer )

    def test_united_group_4(self):

        # preparation
        atom_info = dict(
            names = ['H1','H2','H3','N','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','N','H','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','N','H','CA','HA','CB','HB1','HB2','HB3',
                     'C','O','OXT'],

            ids   = range(1,33+1),

            elems = ['H','H','H','N','C','H','C','H','H','H',
                     'C','O','N','H','C','H','C','H','H','H',
                     'C','O','N','H','C','H','C','H','H','H',
                     'C','O','O'],
        )

        # answer value
        answer = [ ('00004_N', [1, 2, 3, 4]),
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
                   ('00033_OXT', [33])
        ]

        # test
        group_atoms_pairs = group.gen_united_group(atom_info)

        # check
        eq_( list(group_atoms_pairs), answer )

