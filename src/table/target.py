from __future__ import print_function

import numpy

import interact_table as it

def parse_target_atoms_line(target_atoms_line, natom):

    atoms = []

    for mask in target_atoms_line:
        if '-' in mask:
            atoms.extend( parse_mask_with_range(mask, atoms, natom) )

        else:
            atoms.append( int(mask) )
            # print('The format of target_atoms keyword is invalid : ', mask)
            # raise 'Error'

    return atoms

def parse_mask_with_range(mask, atoms, natom):

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
        iatm_end = natom

    else:
        pass

    return range(iatm_beg, iatm_end+1)

def make_iatm_to_itar(target_atoms, natom):
    """Make iatm => itarget dictionary."""
    iatm_to_itar = numpy.zeros( [natom], numpy.int)

    for itar_1, iatm in enumerate(target_atoms):
        iatm_to_itar[iatm-1] = itar_1 + 1

    return iatm_to_itar

def make_interaction_table_current(base_table, target_atoms, natom):
    """Make the interaction table selected by target atoms for current."""

    iatm_to_itar = make_iatm_to_itar(target_atoms, natom)
    iatm_min = min(target_atoms)
    iatm_max = max(target_atoms)

    table = it.InteractionTable(base_table=base_table)

    def squeeze_target_atoms_range(iatm, jatm):
        itar = iatm_to_itar[iatm-1]
        jtar = iatm_to_itar[jatm-1]
        return (itar != 0) or (jtar != 0)

    return table.filter_or(iatm_min, iatm_max).filter_fun(
            squeeze_target_atoms_range)

def make_interaction_table_flux(base_table, target_atoms, natom):
    """Make the interaction table selected by target atoms for flux."""

    iatm_to_itar = make_iatm_to_itar(target_atoms, natom)
    iatm_min = min(target_atoms)
    iatm_max = max(target_atoms)

    table = it.InteractionTable(base_table=base_table)

    def squeeze_target_atoms_range(iatm, jatm):
        itar = iatm_to_itar[iatm-1]
        jtar = iatm_to_itar[jatm-1]
        return (itar != 0) and (jtar != 0)

    return table.filter_and(iatm_min, iatm_max).filter_fun(
            squeeze_target_atoms_range)

if __name__ == '__main__':
    natom = 33
    maskline = ['4',  '7-9', '13-18']
    target_atoms = parse_target_atoms_line(maskline, natom)
    print(target_atoms)

    maskline = ['-15']
    target_atoms = parse_target_atoms_line(maskline, natom)
    print(target_atoms)

    maskline = ['-5', '8-']
    target_atoms = parse_target_atoms_line(maskline, natom)
    print(target_atoms)
