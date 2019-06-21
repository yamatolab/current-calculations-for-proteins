from __future__ import print_function

import os, sys
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)

# additional package

# curp package

################################################################################
class TopologyParser:

    """
    >>> tpl = TopologyParser('./test/ala3.tpl.gz')
    >>> tpl.parse()
    >>> tpl.get_molcule_info()

    """

    section_end = 'TPL> END SECTION'
    comment = ';'

    def __init__(self, filename):
        self._filename = filename

    def get_molcule_info(self):
        """Get name dictionary, name, mol_names, and num_molecules."""
        return self.__mol_info

    def get_atom_info(self, molname=None):
        """Get dictionary of ids, names, types, type_ids, ren_names,
        res_ids, masses, vdw_radii, and charges by molname.
        """
        if molname is None:
            molname = self.get_molcule_info()[0]['name']
        return self.__atom_infos.get(molname)

    def get_bond_info(self, molname=None):
        """Get dictionary of two_atoms, force_consts, and length_eqs
        by molname.
        """
        if molname is None:
            molname = self.get_molcule_info()[0]['name']
        return self.__bond_infos.get(molname)

    def get_angle_info(self, molname=None):
        """Get dictionary of three_atoms, force_consts, and theta_eqs
        by molname.
        """
        if molname is None:
            molname = self.get_molcule_info()[0]['name']
        return self.__angle_infos.get(molname)

    def get_torsion_info(self, molname=None):
        """Get dictionary of four_atoms, force_consts, num_torsions, num_freqs,
        initial_phases, and nonbonded_flags by molname.
        """
        if molname is None:
            molname = self.get_molcule_info()[0]['name']
        return self.__torsion_infos.get(molname)

    def get_improper_info(self, molname=None):
        """Get dictionary of four_atoms, force_consts, num_torsions, num_freqs,
        initial_phases, and nonbonded_flags by molname.
        """
        if molname is None:
            molname = self.get_molcule_info()[0]['name']
        return self.__improper_infos.get(molname)

    def get_nonbond_info(self):
        """Get dictionary of atom_types, vdw_radii and vdw_depths."""
        return self.__nonbond_info

    def parse(self):
        tpl_file = open(self._filename, 'r')
        sections = self.split_content(tpl_file)
        tpl_file.close()

        atom_infos     = {}
        bond_infos     = {}
        angle_infos    = {}
        torsion_infos  = {}
        improper_infos = {}

        mol_info = self.parse_mol(sections['mol'])
        self.print_info(mol_info)
        self.__mol_info = mol_info

        for atom_lines in sections['atom']:
            info = self.parse_atom(atom_lines)
            atom_infos[info['molname']] = info
            self.print_info(info)
        self.__atom_infos = atom_infos

        for bond_lines in sections['bond']:
            info = self.parse_bond(bond_lines)
            bond_infos[info['molname']] = info
            self.print_info(info)
        self.__bond_infos = bond_infos

        for angle_lines in sections['angle']:
            info = self.parse_angle(angle_lines)
            angle_infos[info['molname']] = info
            self.print_info(info)
        self.__angle_infos = angle_infos

        for torsion_lines in sections['torsion']:
            info = self.parse_torsion(torsion_lines)
            torsion_infos[info['molname']] = info
            self.print_info(info)
        self.__torsion_infos = torsion_infos

        for improper_lines in sections['improper']:
            info = self.parse_torsion(improper_lines)
            improper_infos[info['molname']] = info
            self.print_info(info)
        self.__improper_infos = improper_infos

        nonbond_info = self.parse_nonbond(sections['nonbond'])
        self.print_info(nonbond_info)
        self.__nonbond_info = nonbond_info

    def split_content(self, tpl_file):

        atom_lines_list     = []
        bond_lines_list     = []
        angle_lines_list    = []
        torsion_lines_list  = []
        improper_lines_list = []

        gen_lines = self.gen_optimized_lines(tpl_file)
        sections = []
        for line in gen_lines:
            if not self.is_section(line): continue

            lines = list( self.gen_section_lines(gen_lines) )

            if line.startswith('TPL> MOLECULES'):
                self.print_section('[MOLECULES]', lines)
                mol_lines = lines

            elif line.startswith('TPL> ATOMS'):
                self.print_section('[ATOMS]', lines)
                atom_lines_list.append(lines)

            elif line.startswith('TPL> BONDS'):
                self.print_section('[BONDS]', lines)
                bond_lines_list.append(lines)

            elif line.startswith('TPL> ANGLES'):
                self.print_section('[ANGLES]', lines)
                angle_lines_list.append(lines)

            elif line.startswith('TPL> TORSIONS'):
                self.print_section('[TORSIONS]', lines)
                torsion_lines_list.append(lines)

            elif line.startswith('TPL> IMPROPER-TORSIONS'):
                self.print_section('[IMPROPER]', lines)
                improper_lines_list.append(lines)

            # elif line.startswith('TPL> FUNCTIONS'):
            #     self.print_section('FUNCTIONS', lines)
            #     function_lines = lines

            elif line.startswith('TPL> NONBONDS'):
                self.print_section('[NONBONDS]', lines)
                nonbond_lines = lines

            else:
                continue

        return dict(mol=mol_lines, 
                atom=atom_lines_list, bond=bond_lines_list,
                angle=angle_lines_list, torsion=torsion_lines_list,
                improper=improper_lines_list,
                # function=function_lines,
                nonbond=nonbond_lines )

    def gen_optimized_lines(self, lines):
        """Generate lines with the end of sections and without the comments."""

        # search first section marker
        for line in lines:
            line = line.split(self.comment)[0].strip()
            if self.is_section(line):
                yield line
                break

        # rest
        for line in lines:
            line = line.split(self.comment)[0].strip()
            if self.is_section(line):
                yield self.section_end
                yield line
            elif line != '':
                yield line
            else:
                continue
        else:
            yield self.section_end

    def is_section(self, line):
        if line == self.section_end:
            return False
        elif line.startswith('TPL>'):
            return True
        else:
            return False

    def gen_section_lines(self, gen):
        for line in gen:
            if line == self.section_end:
                break
            yield line

    def print_section(self, header, lines):
        logger.debug(header)
        for l in lines:
            logger.debug(l)

    def print_info(self, info):
        """Print the content of information."""
        if 'molname' in info:
            logger.debug('')
            logger.debug('*** MOLCULE NAME: ' + info['molname'] + ' ***')
            logger.debug('')

        for key, values in info.items():
            logger.debug(key)
            logger.debug(values)

    def parse_mol(self, lines):

        mol_names = []
        num_molecules = []

        for line in lines:
            cols = line.split()
            mol_names.append( cols[0] )
            num_molecules.append( int(cols[1]) )

        return dict(name=None,
                mol_names=mol_names, num_molecules=num_molecules )

    def parse_atom(self, lines):

        molname = lines[0]

        ids       = []
        names     = []
        types     = []
        type_ids  = []
        res_names = []
        res_ids   = []
        masses    = []
        vdw_radii = []
        charges   = []

        arrow_flag = False

        for line in lines[1:]:
            if arrow_flag:
                if not line.strip().endswith('->'):
                    arrow_flag = False
                continue
                
            cols = line.split()

            names.append(cols[0])
            types.append(cols[1])
            type_ids.append(cols[2])
            res_names.append(cols[3])
            res_ids.append(int(cols[4]))
            masses.append(float(cols[5]))
            vdw_radii.append(float(cols[6]))
            charges.append(float(cols[7]))

            if line.strip().endswith('->'):
                arrow_flag = True

        return dict(molname=molname,
                ids=ids, names=names, types=types, type_ids=type_ids,
                res_names=res_names, res_ids=res_ids, masses=masses,
                vdw_radii=vdw_radii, charges=charges )

    def parse_bond(self, lines):

        molname = lines[0]

        two_atoms    = []
        force_consts = []
        length_eqs   = []

        for line in lines[1:]:
            cols = line.split()

            two_atoms.append( [int(cols[icol]) for icol in range(2)] )
            force_consts.append( float(cols[2]) ) # unit: kcal/(mol*A*A)
            length_eqs.append( float(cols[3]) )    # unit: A

        return dict(molname=molname,
                two_atoms=two_atoms, force_consts=force_consts,
                length_eqs=length_eqs )

    def parse_angle(self, lines):

        molname = lines[0]

        three_atoms     = []
        force_consts  = []
        theta_eqs = []

        for line in lines[1:]:
            cols = line.split()

            three_atoms.append( [int(cols[icol]) for icol in range(3)] )
            force_consts.append( float(cols[3]) ) # kcal/(mol*rad*rda)
            theta_eqs.append( float(cols[4]) ) # degree

        return dict(molname=molname,
                three_atoms=three_atoms, force_consts=force_consts,
                theta_eqs=theta_eqs )

    def parse_torsion(self, lines):

        molname = lines[0]

        four_atoms     = []
        num_torsions   = []
        num_freqs      = []
        force_consts   = []
        initial_phases = []
        nonbonded_flags = []

        for line in lines[1:]:
            cols = line.split()

            four_atoms.append( [int(cols[icol]) for icol in range(4)] )
            force_consts.append( float(cols[4]) )
            num_torsions.append( int(cols[5]) )
            num_freqs.append( int(cols[6]) )
            initial_phases.append( float(cols[7]) )
            nonbonded_flags.append( int(cols[8]) )

        return dict(molname=molname,
                four_atoms=four_atoms, force_consts=force_consts,
                num_torsions=num_torsions, num_freqs=num_freqs,
                initial_phases=initial_phases, 
                nonbonded_flags=nonbonded_flags )

    def parse_nonbond(self, lines):

        atom_types = []
        vdw_radii  = []
        vdw_depths = []

        for line in lines:
            cols = line.split()
            if len(cols) != 7: continue

            atom_types.append( int(cols[0]) )
            vdw_radii.append( float(cols[3]) )
            vdw_depths.append( float(cols[4]) )

        return dict(title=None,
                atom_types=atom_types, vdw_radii=vdw_radii,
                vdw_depths=vdw_depths )


################################################################################

from forcefield.amberbase import ConverterBase
class Format2Amber99Converter(ConverterBase):

    """
    A data in PRESTO version 2 style topology => force field
    """

    def __init__(self, topology):
        ConverterBase.__init__(self, topology)
        self.__mol_info = []
        self.__residue_info = {}
        self.__nonbond_interact_table = None

        self.__bond_info     = {}
        self.__angle_info    = {}
        self.__torsion_info  = {}
        self.__improper_info = {}
        self.__coulomb_info  = {}
        self.__vdw_info      = {}

        self.__14bonded_pairs = None

    def convert(self):
        # bonded interaction
        self._convert_residue()
        self._convert_bond()
        self._convert_angle()
        self._convert_torsion()
        self._convert_improper()
        self._convert_coulomb()
        self._convert_vdw()

    def get_mol_info(self):
        if not self.__mol_info:
            self._store_mol_info()
        return self.__mol_info

    def get_residue_info(self):
        return self.__residue_info

    def get_bond_info(self):
        return self.__bond_info

    def get_angle_info(self):
        return self.__angle_info

    def get_torsion_info(self):
        return self.__torsion_info

    def get_improper_info(self):
        return self.__improper_info

    def get_coulomb_info(self):
        return self.__coulomb_info

    def get_vdw_info(self):
        return self.__vdw_info

    def _store_mol_info(self):
        mol_info = self.get_topology().get_molcule_info()
        for name, nmol in zip(mol_info['mol_names'], mol_info['num_molecules']):
            natom = len(self.get_topology().get_atom_info(name)['names'])
            new_info = dict(name=name, nmol=nmol, natom=natom)
            self.__mol_info.append(new_info)

    def get_natom(self):
        natom = 0
        for info in self.get_mol_info():
            natom += info['nmol'] * info['natom']
        return natom

    def convert_atoms(self, atoms_info):

        self.topology.store_params('coulomb', charges=charges)

    def _convert_residue(self):

        ids       = []
        names     = []
        nres_last = 0

        # loop for the number of molcule kinds
        for minfo in self.get_mol_info():

            atom_info = self.get_topology().get_atom_info(minfo['name'])
            if atom_info:

                # loop for the number of molecules
                for imol in range( minfo['nmol'] ):
                    names.extend(atom_info['res_names'])

                    ids.extend( resid + nres_last
                            for resid in atom_info['res_ids'] )
                    nres_last = ids[-1]

        self.__residue_info = dict(ids=ids, names=names)

    def _convert_bond(self):

        two_atoms    = []
        force_consts = []
        length_eqs   = []
        cur_natom = 0

        # loop for the number of molcule kinds
        for minfo in self.get_mol_info():

            bnd_info = self.get_topology().get_bond_info(minfo['name'])
            if bnd_info:

                # loop for the number of molecules
                for imol in range( minfo['nmol'] ):

                    for iatoms in bnd_info['two_atoms']:
                        two_atoms.append(
                                [iatoms[i]+cur_natom for i in range(2)] )

                    force_consts.extend(bnd_info['force_consts'])
                    length_eqs.extend(bnd_info['length_eqs'])

                    cur_natom += minfo['natom']

            else:
                cur_natom += minfo['natom'] * minfo['nmol']

        self.__bond_info = dict(
                two_atoms=two_atoms,
                force_consts=force_consts,
                length_eqs=length_eqs )

    def _convert_angle(self):

        three_atoms  = []
        force_consts = []
        theta_eqs    = []
        cur_natom    = 0

        # loop for the number of molcule kinds
        for minfo in self.get_mol_info():

            ang_info = self.get_topology().get_angle_info(minfo['name'])
            if ang_info:

                # loop for the number of molecules
                for imol in range( minfo['nmol'] ):

                    for iatoms in ang_info['three_atoms']:
                        three_atoms.append(
                                [iatoms[i]+cur_natom for i in range(3)] )

                    force_consts.extend(ang_info['force_consts'])
                    theta_eqs.extend(ang_info['theta_eqs'])
                    cur_natom += minfo['natom']

            else:
                cur_natom += minfo['natom'] * minfo['nmol']

        self.__angle_info = dict(
                three_atoms=three_atoms,
                force_consts=force_consts,
                theta_eqs=theta_eqs )

    def _convert_torsion(self):
        self.__torsion_info  = self.__convert_torsion(
                self.get_topology().get_torsion_info)

    def _convert_improper(self):
        self.__improper_info = self.__convert_torsion(
                self.get_topology().get_improper_info)

    def __convert_torsion(self, get_info_fun):

        four_atoms     = []
        num_torsions   = []
        num_freqs      = []
        force_consts   = []
        initial_phases = []
        cur_natom = 0

        # loop for the number of molcule kinds
        for minfo in self.get_mol_info():

            tor_info = get_info_fun(minfo['name'])
            if tor_info:

                # loop for the number of molecules
                for imol in range( minfo['nmol'] ):

                    for iatoms in tor_info['four_atoms']:
                        four_atoms.append(
                                [iatoms[i]+cur_natom for i in range(4)] )

                    num_torsions.extend(tor_info['num_torsions'])
                    num_freqs.extend(tor_info['num_freqs'])
                    force_consts.extend(tor_info['force_consts'])
                    initial_phases.extend(tor_info['initial_phases'])

                    cur_natom += minfo['natom']

            else:
                cur_natom += minfo['natom'] * minfo['nmol']

        all_torsion_info = dict(
                four_atoms=four_atoms,
                num_torsions=num_torsions,
                num_freqs=num_freqs,
                force_consts=force_consts,
                initial_phases=initial_phases )

        return all_torsion_info

    def _convert_coulomb(self):

        charges   = []
        cur_natom = 0

        # loop for the number of molcule kinds
        for minfo in self.get_mol_info():

            atom_info = self.get_topology().get_atom_info(minfo['name'])
            if atom_info:

                # loop for the number of molecules
                for imol in range( minfo['nmol'] ):
                    charges.extend(atom_info['charges'])
                    cur_natom += minfo['natom']

            else:
                cur_natom += minfo['natom'] * minfo['nmol']

        self.__coulomb_info = dict(charges=charges)

    def _convert_vdw(self):

        cur_natom = 0
        atom_types = []

        # atom_types
        # loop for the number of molcule kinds
        for minfo in self.get_mol_info():

            atom_info = self.get_topology().get_atom_info(minfo['name'])
            if atom_info:

                # loop for the number of molecules
                for imol in range( minfo['nmol'] ):
                    atom_types.extend(atom_info['type_ids'])
                    cur_natom += minfo['natom']

            else:
                cur_natom += minfo['natom'] * minfo['nmol']

        # c6s and c12s
        """Get dictionary of atom_types, vdw_radii and vdw_depths."""
        nonbond_info = self.get_topology().get_nonbond_info()

        types_radii_depths = zip(
                nonbond_info['atom_types'],
                nonbond_info['vdw_radii'],
                nonbond_info['vdw_depths'] )

        ntype = len(nonbond_info['atom_types'])
        c6s =  numpy.zeros((ntype,ntype))
        c12s = numpy.zeros((ntype,ntype))
        
        import math
        for itype, rad_i, dep_i in types_radii_depths:
            for jtype, rad_j, dep_j in types_radii_depths:
                dep_ij = math.sqrt(dep_i * dep_j)
                rad_ij = rad_j + rad_j
                c6s[itype-1, jtype-1]  = 2.0 * dep_ij * rad_ij**6
                c12s[itype-1, jtype-1] = dep_ij * rad_ij**12

        self.__vdw_info = dict(atom_types=atom_types, c6s=c6s, c12s=c12s)

    def _make_14bonded_pairs(self):

        cur_natom = 0
        two_atoms = []

        # loop for the number of molcule kinds
        for minfo in self.get_mol_info():

            tor_info = self.get_topology().get_torsion_info(minfo['name'])
            if tor_info:

                # loop for the number of molecules
                for imol in range( minfo['nmol'] ):

                    atoms_flags = zip(
                            tor_info['nonbonded_flags'],
                            tor_info['four_atoms'] )

                    ####################### flag ###############################
                    two_atoms.extend(
                        [ (iatoms[0]+cur_natom, iatoms[3]+cur_natom) 
                            for flag, iatoms in atoms_flags if flag==1 ] )
                    ####################### flag ###############################

                    cur_natom += minfo['natom']

            else:
                cur_natom += minfo['natom'] * minfo['nmol']

        self.__14bonded_pairs = two_atoms

    def get_14bonded_pairs(self):
        if self.__14bonded_pairs is None:
            self._make_14bonded_pairs()
        return self.__14bonded_pairs
    get_bonded14_pairs = get_14bonded_pairs


################################################################################
class Format2OPLSConverter(Format2Amber99Converter):

    def _convert_vdw(self):

        cur_natom = 0
        interact_table = None
        atom_types = []

        # atom_types
        # loop for the number of molcule kinds
        for minfo in self.get_mol_info():

            atom_info = self.get_topology().get_atom_info(minfo['name'])
            if atom_info:

                # loop for the number of molecules
                for imol in range( minfo['nmol'] ):
                    atom_types.extend(atom_info['type_ids'])
                    cur_natom += minfo['natom']

            else:
                cur_natom += minfo['natom'] * minfo['nmol']

        # c6s and c12s
        """Get dictionary of atom_types, vdw_radii and vdw_depths."""
        nonbond_info = self.get_topology().get_nonbond_info()

        types_radii_depths = zip(
                nonbond_info['atom_types'],
                nonbond_info['vdw_radii'],
                nonbond_info['vdw_depths'] )

        ntype = len(nonbond_info['atom_types'])
        c6s =  numpy.zeros((ntype,ntype))
        c12s = numpy.zeros((ntype,ntype))
        
        import math
        for itype, rad_i, dep_i in types_radii_depths:
            for jtype, rad_j, dep_j in types_radii_depths:
                dep_ij = math.sqrt(dep_i * dep_j)
                rad_ij = rad_j * rad_j
                c6s[itype-1, jtype-1]  = 4.0 * dep_ij * rad_ij**3
                c12s[itype-1, jtype-1] = 4.0 * dep_ij * rad_ij**6

        all_vdw_info = dict( atom_types=atom_types, c6s=c6s, c12s=c12s )

        logger.debug('*** vdw ***')
        for key, values in all_vdw_info.items():
            logger.debug('[{}]'.format(key))
            logger.debug(values)

        return all_vdw_info


if __name__ == '__main__':
     logger.log_level = logger.INFO
     tpl = TopologyParser('./test/pyp.tpl')
     tpl.parse()

     converter = Format2Amber99Converter(tpl)
     converter.convert()
     converter.print_residue()
     converter.print_bond()
     converter.print_angle()
     converter.print_torsion()
     converter.print_improper()
     converter.print_coulomb()
     converter.print_vdw()
     converter.print_bonded_pairs()
     converter.print_14bonded_pairs()
