from __future__ import print_function

import os, sys
from collections import OrderedDict as odict
import numpy

topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)
import exception
import clog as logger

class SectionNotFoundError(exception.CurpException): pass
class ParserError(exception.CurpException): pass


def get_group_iatoms_pairs(setting, target_atoms, res_info, atom_info):
    """Get group atoms table with reference to setting."""

    method = setting.curp.group_method

    natom = len( atom_info['elems'] )
    iatm_to_itars = get_iatm_to_itars(target_atoms, natom)

    if method == 'none':
        gname_iatoms_pairs = []

    elif method == 'residue':
        gname_iatoms_pairs = list( gen_group_pairs_with_target(
            gen_residue_group(res_info, atom_info), iatm_to_itars ) )

    elif method == 'united':
        gname_iatoms_pairs = list( gen_group_pairs_with_target(
            gen_united_group(atom_info), iatm_to_itars ) )

    elif method == 'file':
        grp_parser = GroupParser(natom, setting.input.group_file[0])
        gname_iatoms_pairs = list( gen_group_pairs_with_target(
            grp_parser.gen_group_atoms_pairs(), iatm_to_itars ) )

    else:
        raise Exception

    logger.info_title('Group information')
    for gname, iatoms in gname_iatoms_pairs:
        logger.info('[',gname,']')
        logger.info( ' '.join(str(iatm) for iatm in iatoms) )

    return gname_iatoms_pairs


import ini_parser as ini
class GroupParser(ini.IniParser):

    def __init__(self, natom, *args, **kwds):
        self.__natom = natom
        ini.IniParser.__init__(self, *args, **kwds)

    def gen_group_atoms_pairs(self):
        return ini.IniParser.gen_secname_cols_pair(self)

    def _gen_col(self, lines):
        atoms = []
        for line in lines:
            if '=' in line: continue
            for col in line.split():

                if '-' in col:
                    atoms.extend( self._parse_mask_with_range(col, atoms) )

                else:
                    atoms.append( int(col) )

        return atoms

    def _parse_mask_with_range(self, mask, atoms):

        cols = mask.split('-')

        if cols[0] != '' and cols[1] != '':
            iatm_beg, iatm_end = int(cols[0]), int(cols[1])

        elif cols[0] == '' and cols[1] != '': # -45
            iatm_end = int(cols[1])

            if len(atoms) == 0:
                iatm_beg = 1
            else:
                iatm_beg = atoms[-1] + 1

        elif cols[0] != '' and cols[1] == '': # 45-
            iatm_beg = int(cols[0])
            iatm_end = self.__natom

        else:
            pass

        return range(iatm_beg, iatm_end+1)

# for topology and bonded
def squeeze_calculating_bonded(groups, table_name, atoms_table, **other_table):
    """Squeeze the bonded calculation items in crude and fast."""

    # get minimum and maximum
    gmin = min( min(iatoms) for iatoms in groups.values() )
    gmax = max( max(iatoms) for iatoms in groups.values() )

    # make index list
    squeezed_indexes = []
    for index, iatoms in enumerate(atoms_table):
        for iatm in iatoms:
            if gmin <= iatm <= gmax:
                squeezed_indexes.append(index)
                break
            else:
                continue

    # reconstruct data
    results = {}
    results[table_name] = [ atoms_table[index] for index in squeezed_indexes ]
    for key, values in other_table.items():
        results[key] = [ values[index] for index in squeezed_indexes ]

    return results


def squeeze_calculating_bonded2(groups, table_name, atoms_table, **other_table):
    """Squeeze the bonded calculation items in precise and fast."""

    # get flatten group list
    group_atoms = []
    for atoms in groups.values():
        group_atoms.extend(atoms)

    # make index list
    squeezed_indexes = []
    for index, iatoms in enumerate(atoms_table):
        for iatm in iatoms:
            if iatm in group_atoms:
                squeezed_indexes.append(index)
                break
            else:
                continue

    # reconstruct data
    results = {}
    results[table_name] = [ atoms_table[index] for index in squeezed_indexes ]
    for key, values in other_table.items():
        results[key] = [ values[index] for index in squeezed_indexes ]

    return results


def gen_residue_group(res_info, atom_info, name_fmt='', rid_first=1):

    name_fmt = name_fmt if name_fmt else '{rid:05}_{rname}'
    res_ids = res_info['ids']
    res_names = res_info['names']
    natom = len(atom_info['names'])

    rid_last = rid_first

    iatoms = []
    for iatm_1, rid in zip(range(natom), res_ids):

        iatm = iatm_1 + 1

        if rid != rid_last:
            rname = res_names[iatm_1-1]
            gname = name_fmt.format(rname=rname, rid=rid-1)
            yield (gname, iatoms)

            rid_last = rid
            iatoms = []

        iatoms.append( iatm )

    else:
        rname = res_names[-1]
        rid = res_ids[-1]
        gname = name_fmt.format(rname=rname, rid=rid)
        yield (gname, iatoms)
    

def gen_united_group(atom_info, name_fmt=''):

    name_fmt = name_fmt if name_fmt else '{id:05}_{name}'

    # atom_ids = atom_info['ids']
    atom_elems = atom_info['elems']
    atom_names = atom_info['names']
    natom = len(atom_elems)

    # first step
    iatoms = []
    i_beg = 0
    for iatm_1, (elem, name) in enumerate(zip(atom_elems, atom_names)):
        iatm = iatm_1 + 1

        i_beg += 1
        if elem != 'H':
            name_last = name
            iatm_last = iatm
            iatoms.append( iatm )
            break

        iatoms.append( iatm )

    for iatm_1, elem, name in zip(
            range(i_beg, natom+1), atom_elems[i_beg:], atom_names[i_beg:] ):
        iatm = iatm_1 + 1

        if elem != 'H':
            gname = name_fmt.format(name=name_last, id=iatm_last)
            yield (gname, iatoms)
            iatoms = [iatm]
            name_last = name
            iatm_last = iatm

        else:
            iatoms.append( iatm )

    else:
        gname = name_fmt.format(name=name_last, id=iatm_last)
        yield (gname, iatoms)


def gen_group_pairs_with_target(gen_gname_iatoms_pairs, iatm_to_itars):

    for gname, iatoms in gen_gname_iatoms_pairs:

        iatoms_new = []
        for iatm in iatoms:
            if iatm_to_itars[iatm-1] == 0: continue
            iatoms_new.append( iatm )

        if len(iatoms_new) != 0:
            yield (gname, iatoms_new)


def get_iatm_to_itars(target_atoms, natom):
    iatm_to_itars = numpy.zeros( [natom], numpy.int)
    for itar_1, iatm in enumerate(target_atoms):
        iatm_to_itars[iatm-1] = itar_1 + 1
    return iatm_to_itars

