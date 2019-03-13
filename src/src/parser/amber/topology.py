from __future__ import print_function

import os, sys
import gzip
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)

# additional package
import numpy

# curp package
import table.interact_table as it
import clog as logger

################################################################################
class TopologyParser:

    """
    >>> tpl = TopologyParser('./test/ala3.prmtop')
    >>> tpl.parse()
    >>> tpl.get_molcule_info()

    """

    section_end = '%END FLAG'
    comment = ';'

    def __init__(self, filename):
        self._filename = filename
        self.__fname_to_info = {}

    def get_param_infos(self, flag_name):
        return self.__fname_to_info[flag_name]

    def parse(self):
        if self._filename.endswith('.gz'):
            tpl_file = gzip.open(self._filename, 'rb')
        else:
            tpl_file = open(self._filename, 'r')
        fname_to_lines = self.split_content(tpl_file)
        tpl_file.close()

        fname_to_parsers = {
            'ATOM_NAME'                  : self._parse_name(20, 4)    , 
            'CHARGE'                     : self._parse_generic(float) , 
            'MASS'                       : self._parse_generic(float) , 
            'ATOM_TYPE_INDEX'            : self._parse_generic(int)   , 
            'RESIDUE_LABEL'              : self._parse_name(20, 4)    , 
            'RESIDUE_POINTER'            : self._parse_generic(int)   , 
            'BOND_FORCE_CONSTANT'        : self._parse_generic(float) , # kcal/mol
            'BOND_EQUIL_VALUE'           : self._parse_generic(float) , # A
            'BONDS_INC_HYDROGEN'         : self._parse_generic(int)   , 
            'BONDS_WITHOUT_HYDROGEN'     : self._parse_generic(int)   , 
            'ANGLE_FORCE_CONSTANT'       : self._parse_generic(float) , # kcal/mol A^2
            'ANGLE_EQUIL_VALUE'          : self._parse_generic(float) , # radian
            'ANGLES_INC_HYDROGEN'        : self._parse_generic(int)   , 
            'ANGLES_WITHOUT_HYDROGEN'    : self._parse_generic(int)   , 
            'DIHEDRAL_FORCE_CONSTANT'    : self._parse_generic(float) , # kcal/mol
            'DIHEDRAL_PERIODICITY'       : self._parse_generic(float) , 
            'DIHEDRAL_PHASE'             : self._parse_generic(float) , 
            'DIHEDRALS_INC_HYDROGEN'     : self._parse_generic(int)   , 
            'DIHEDRALS_WITHOUT_HYDROGEN' : self._parse_generic(int)   , 
            'NONBONDED_PARM_INDEX'       : self._parse_generic(int)   , 
            'LENNARD_JONES_ACOEF'        : self._parse_generic(float) , 
            'LENNARD_JONES_BCOEF'        : self._parse_generic(float) , 
            'HBOND_ACOEF'                : self._parse_generic(float) , 
            'HBOND_BCOEF'                : self._parse_generic(float) , 
            'AMBER_ATOM_TYPE'            : self._parse_name(20, 4)    , 
            'BOX_DIMENSIONS'             : self._parse_generic(float) , 
        }
        for fname, parser in fname_to_parsers.items():
            if fname in fname_to_lines:
                lines = fname_to_lines[fname]
                self.__fname_to_info[fname] = list(parser(lines))
            else:
                self.__fname_to_info[fname] = None

        for fname, info in self.__fname_to_info.items():
            self.print_info(fname, info)

    def split_content(self, tpl_file):

        fname_to_lines = {}

        gen_lines = self.gen_optimized_lines(tpl_file)
        sections = []
        for line in gen_lines:
            if not self.is_section(line): continue

            flagname = line.split()[1]
            lines = list( self.gen_section_lines(gen_lines) )
            fname_to_lines[flagname] = lines 

        return fname_to_lines

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
            elif line.startswith('%FORMAT'):
                continue
            elif line != '':
                yield line
            else:
                continue
        else:
            yield self.section_end

    def is_section(self, line):
        if line == self.section_end:
            return False
        elif line.startswith('%FLAG'):
            return True
        else:
            return False

    def gen_section_lines(self, gen):
        for line in gen:
            if line == self.section_end:
                break
            yield line

    def print_info(self, header, data):
        logger.debug('['+header+']')
        logger.debug(data)

    def _parse_generic(self, convtype):
        def parse(lines):
            for line in lines:
                cols = line.split()
                for col in cols:
                    yield convtype(col)
        return parse

    def _parse_name(self, ncol=20, width=4):
        def parse(lines):
            for line in lines:
                for icol in range(ncol):
                    name = line[icol*width:(icol+1)*width].strip()
                    if name == '': continue
                    yield name
        return parse


################################################################################


# definition of vdw radius parameters.
type_to_vdws = {
        'H'  : 0.6000 , 
        # 'HO' : 0.0000 , 
        'HO' : 0.6000 ,  # Artifact
        'HS' : 0.6000 , 
        'HC' : 1.4870 , 
        'H1' : 1.3870 , 
        'H2' : 1.2870 , 
        'H3' : 1.1870 , 
        'HP' : 1.1000 , 
        'HA' : 1.4590 , 
        'H4' : 1.4090 , 
        'H5' : 1.3590 , 
        # 'HW' : 0.0000 , 
        'HW' : 0.6000 , # Artifact
        'HZ' : 1.4590 , 
        'O'  : 1.6612 , 
        'O2' : 1.6612 , 
        'OW' : 1.7683 , 
        'OH' : 1.7210 , 
        'OS' : 1.6837 , 
        'C*' : 1.9080 , 
        'CT' : 1.9080 , 
        'CX' : 1.9080 , 
        'C'  : 1.9080 , 
        'N'  : 1.8240 , 
        'N3' : 1.8240 , 
        'NY' : 1.8240 , 
        'S'  : 2.0000 , 
        'SH' : 2.0000 , 
        'P'  : 2.1000 , 
        'IM' : 2.47   , 
        'Li' : 1.1370 , 
        'IP' : 1.8680 , 
        'Na' : 1.8680 , 
        'Na+': 1.8680 , 
        'K'  : 2.6580 , 
        'Rb' : 2.9560 , 
        'Cs' : 3.3950 , 
        'MG' : 0.7926 , 
        'C0' : 1.7131 , 
        'Zn' : 1.10   , 
        'F'  : 1.75   , 
        'Cl' : 1.948  , 
        'Cl-': 1.948  , 
        'Br' : 2.22   , 
        'I'  : 2.35   , 
        'IB' : 5.0    , 
        'LP' : 0.00   , 
  }

# for gaff
gaff_type_to_vdws = {
        'h1' : 1.3870 , 
        'h2' : 1.2870 , 
        'h3' : 1.1870 , 
        'h4' : 1.4090 , 
        'h5' : 1.3590 , 
        'ha' : 1.4590 , 
        'hc' : 1.4870 , 
        'hn' : 0.6000 , 
        'ho' : 0.0000 , 
        'hp' : 0.6000 , 
        'hs' : 0.6000 , 
        'hw' : 0.0000 , 
        'hx' : 1.1000 , 
        'o'  : 1.6612 , 
        'oh' : 1.7210 , 
        'os' : 1.6837 , 
        'ow' : 1.7683 , 
        'c'  : 1.9080 , 
        'c1' : 1.9080 , 
        'c2' : 1.9080 , 
        'c3' : 1.9080 , 
        'ca' : 1.9080 , 
        'cc' : 1.9080 , 
        'cd' : 1.9080 , 
        'ce' : 1.9080 , 
        'cf' : 1.9080 , 
        'cg' : 1.9080 , 
        'ch' : 1.9080 , 
        'cp' : 1.9080 , 
        'cq' : 1.9080 , 
        'cu' : 1.9080 , 
        'cv' : 1.9080 , 
        'cx' : 1.9080 , 
        'cy' : 1.9080 , 
        'cz' : 1.9080 , 
        'n'  : 1.8240 , 
        'n1' : 1.8240 , 
        'n2' : 1.8240 , 
        'n3' : 1.8240 , 
        'n4' : 1.8240 , 
        'na' : 1.8240 , 
        'nb' : 1.8240 , 
        'nc' : 1.8240 , 
        'nd' : 1.8240 , 
        'ne' : 1.8240 , 
        'nf' : 1.8240 , 
        'nh' : 1.8240 , 
        'no' : 1.8240 , 
        's'  : 2.0000 , 
        's2' : 2.0000 , 
        's4' : 2.0000 , 
        's6' : 2.0000 , 
        'sx' : 2.0000 , 
        'sy' : 2.0000 , 
        'sh' : 2.0000 , 
        'ss' : 2.0000 , 
        'p2' : 2.1000 , 
        'p3' : 2.1000 , 
        'p4' : 2.1000 , 
        'p5' : 2.1000 , 
        'pb' : 2.1000 , 
        'pc' : 2.1000 , 
        'pd' : 2.1000 , 
        'pe' : 2.1000 , 
        'pf' : 2.1000 , 
        'px' : 2.1000 , 
        'py' : 2.1000 , 
        'f'  : 1.75   , 
        'cl' : 1.948  , 
        'br' : 2.02   , 
        'i'  : 2.15   , 
  }

# for lipid11
lipid_type_to_wdws = {
        'cA' : 1.9080 , 
        'cB' : 1.9080 , 
        'cC' : 1.9080 , 
        'oC' : 1.6612 , 
        'oS' : 1.6837 , 
        'oH' : 1.7210 , 
        'oT' : 1.6837 , 
        'oP' : 1.6612 , 
        'oO' : 1.6612 , 
        'nA' : 1.8240 , 
        'pA' : 2.1000 , 
        'hA' : 1.4870 , 
        'hE' : 1.3870 , 
        'hX' : 1.1000 , 
        'hB' : 1.4590 , 
        'hN' : 0.6000 , 
        'hO' : 0.0000 , 
        'cR' : 1.9080 , 
        'cP' : 1.9080 , 
        'oR' : 1.7210 , 
        'hR' : 0.0000 , 
        'hS' : 1.3870 , 
}

type_to_vdws.update(gaff_type_to_vdws)
type_to_vdws.update(lipid_type_to_wdws)

c_replace_list = ['CA','CB','CC','CD','CK','CM','CN','CQ','CR','CV','CW',
                  'CY','CZ','CP','CS',
                  'CO']
n_replace_list = ['NA','N2','N*','NC','NB','NT','NY']

from forcefield.amberbase import ConverterBase
class Format2AmberBaseConverter(ConverterBase):

    """
    This class converts the Amber style parameter set into the CURP style data.
    """

    def __init__(self, topology, use_atomtype=True):
        ConverterBase.__init__(self, topology)
        self.__mol_info = []
        self.__residue_info = {}
        self.__nonbond_interact_table = None

        self.__atom_info     = {}
        self.__bond_info     = {}
        self.__angle_info    = {}
        self.__torsion_info  = {}
        self.__improper_info = {}
        self.__coulomb_info  = {}
        self.__vdw_info      = {}

        self.__use_atomtype = use_atomtype

    def convert(self, format_setting=None):
        # bonded interaction
        self._convert_residue()
        self._convert_atom()
        self._convert_bond()
        self._convert_angle()
        self._convert_torsion()
        self._convert_improper()
        self._convert_coulomb()
        self._convert_vdw()

        if format_setting is None:
            return

        if format_setting.dump_parameters:
            self.print_atom()
            self.print_residue()
            self.print_bond()
            self.print_angle()
            self.print_torsion()
            self.print_improper()
            self.print_coulomb()
            self.print_vdw()
            self.print_bonded14_pairs()
            self.print_bonded_pairs()

    def get_mol_info(self):
        if not self.__mol_info:
            self._store_mol_info()
        return self.__mol_info

    def get_residue_info(self):
        return self.__residue_info

    def get_atom_info(self):
        return self.__atom_info

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

    def get_pbc_info(self):
        tpl = self.get_topology()
        box_info = tpl.get_param_infos('BOX_DIMENSIONS')
        return box_info

    def _store_mol_info(self):
        mol_info = self.get_topology().get_molcule_info()
        for name, nmol in zip(mol_info['mol_names'], mol_info['num_molecules']):
            natom = len(self.get_topology().get_atom_info(name)['names'])
            new_info = dict(name=name, nmol=nmol, natom=natom)
            self.__mol_info.append(new_info)

    def get_natom(self):
        tpl = self.get_topology()
        atom_names = tpl.get_param_infos('ATOM_NAME')
        return len(atom_names)

    def _convert_residue(self):

        ids       = []
        names     = []
        nres_last = 0
        natom = self.get_natom()

        tpl = self.get_topology()
        amb_res_ptrs = tpl.get_param_infos('RESIDUE_POINTER')
        amb_resnames = tpl.get_param_infos('RESIDUE_LABEL')

        # append the end pointer for each residue
        def gen_detailed_pointer(res_ptrs, natom):
            nres = len(amb_res_ptrs)
            for ires_1 in range(nres):
                iatm_beg = res_ptrs[ires_1]
                if ires_1+1 >= nres:
                    iatm_end = natom
                else:
                    iatm_end = res_ptrs[ires_1+1] - 1
                yield iatm_beg, iatm_end

        name_dptrs = zip(amb_resnames, gen_detailed_pointer(amb_res_ptrs,natom))
        for ires_1, (rname, res_dptr) in enumerate(name_dptrs):
            iatm_beg, iatm_end = res_dptr
            for iatm in range(iatm_beg, iatm_end+1):
                ids.append( ires_1 + 1 )
                names.append( rname )

        self.__residue_info = dict(ids=ids, names=names)

    def _convert_atom(self):
        tpl = self.get_topology()
        names       = tpl.get_param_infos('ATOM_NAME')
        masses      = tpl.get_param_infos('MASS')
        atom_types  = tpl.get_param_infos('AMBER_ATOM_TYPE')

        from math import pi

        elems = []
        for aname in names:
            if aname[0].isdigit():
                elems.append( aname[1] )
            else:
                elems.append( aname[0] )

        if self.__use_atomtype:
            vdw_radii = []
            for i_1, tname in enumerate(atom_types):
                radius = type_to_vdws.get(tname)
                if radius is not None:
                    # print('ok', tname)
                    vdw_radii.append( radius )
                else:

                    if tname in c_replace_list:
                        # print('C* ok')
                        vdw_radii.append( type_to_vdws['C*'] )

                    elif tname in n_replace_list:
                        vdw_radii.append( type_to_vdws['N'] )

                    elif tname[0].isdigit() and tname[1] == 'C':
                        vdw_radii.append( type_to_vdws['C*'] )

                    elif tname[1].isdigit() and tname[0] == 'C':
                        vdw_radii.append( type_to_vdws['C*'] )

                    elif tname[1].isdigit() and tname[0] == 'c':
                        vdw_radii.append( type_to_vdws['C*'] )

                    else:
                        # print('ng', tname)
                        type_to_vdws[tname]

        else:
            vdw_radii = []

        ids = range(1, self.get_natom()+1)
        self.__atom_info = dict(ids=ids, names=names, masses=masses, 
                elems=elems, vdw_radii=vdw_radii)

    def _convert_bond(self):

        tpl = self.get_topology()
        amb_consts    = tpl.get_param_infos('BOND_FORCE_CONSTANT')
        amb_eqs       = tpl.get_param_infos('BOND_EQUIL_VALUE')
        amb_atoms_hyd = tpl.get_param_infos('BONDS_INC_HYDROGEN')
        amb_atoms_woh = tpl.get_param_infos('BONDS_WITHOUT_HYDROGEN')

        two_atoms    = []
        force_consts = []
        length_eqs   = []
        ncol = 3

        # for bonds including hydrogen
        for ibnd in range( len(amb_atoms_hyd)/3 ):
            iatoms  = [ amb_atoms_hyd[ncol*ibnd+i]/3 + 1 for i in range(ncol-1) ]
            index = amb_atoms_hyd[ncol*(ibnd+1) - 1]

            two_atoms.append( iatoms )
            force_consts.append( amb_consts[index-1] )
            length_eqs.append( amb_eqs[index-1] )

        # for bonds without hydrogen
        for ibnd in range( len(amb_atoms_woh)/3 ):
            iatoms = [ amb_atoms_woh[ncol*ibnd+i]/3 + 1 for i in range(ncol-1) ]
            index = amb_atoms_woh[ncol*(ibnd+1) - 1]

            two_atoms.append( iatoms )
            force_consts.append( amb_consts[index-1] )
            length_eqs.append( amb_eqs[index-1] )

        self.__bond_info = dict(
                two_atoms=two_atoms,
                force_consts=force_consts,
                length_eqs=length_eqs )

    def _convert_angle(self):

        tpl = self.get_topology()
        amb_consts    = tpl.get_param_infos('ANGLE_FORCE_CONSTANT')
        amb_eqs       = tpl.get_param_infos('ANGLE_EQUIL_VALUE')
        amb_atoms_hyd = tpl.get_param_infos('ANGLES_INC_HYDROGEN')
        amb_atoms_woh = tpl.get_param_infos('ANGLES_WITHOUT_HYDROGEN')

        three_atoms  = []
        force_consts = []
        theta_eqs    = []
        ncol = 4

        # for angles including hydrogen
        for iang in range( len(amb_atoms_hyd)/ncol ):
            iatoms = [ amb_atoms_hyd[ncol*iang+i]/3 + 1 for i in range(ncol-1) ]
            index = amb_atoms_hyd[ncol*(iang+1) - 1]

            three_atoms.append( iatoms )
            force_consts.append( amb_consts[index-1] )
            theta_eqs.append( amb_eqs[index-1] )

        # for angles without hydrogen
        for iang in range( len(amb_atoms_woh)/4 ):
            iatoms = [ amb_atoms_woh[ncol*iang+i]/3 + 1 for i in range(ncol-1) ]
            index = amb_atoms_woh[ncol*(iang+1) - 1]

            three_atoms.append( iatoms )
            force_consts.append( amb_consts[index-1] )
            theta_eqs.append( amb_eqs[index-1] )

        self.__angle_info = dict(
                three_atoms=three_atoms,
                force_consts=force_consts,
                theta_eqs=theta_eqs )

    def _convert_torsion(self):
        self.__torsion_info = self.__convert_torsion(lambda i: i >= 0)

    def _convert_improper(self):
        self.__improper_info = self.__convert_torsion(lambda i: i < 0)

    def __convert_torsion(self, if_fun):

        tpl = self.get_topology()
        amb_consts    = tpl.get_param_infos('DIHEDRAL_FORCE_CONSTANT')
        amb_freqs     = tpl.get_param_infos('DIHEDRAL_PERIODICITY')
        amb_phases    = tpl.get_param_infos('DIHEDRAL_PHASE')
        amb_atoms_hyd = tpl.get_param_infos('DIHEDRALS_INC_HYDROGEN')
        amb_atoms_woh = tpl.get_param_infos('DIHEDRALS_WITHOUT_HYDROGEN')

        four_atoms     = []
        num_torsions   = []
        num_freqs      = []
        force_consts   = []
        initial_phases = []
        ncol = 5

        # for torsion including hydrogen
        for itor in range( len(amb_atoms_hyd)/ncol ):
            if not if_fun(amb_atoms_hyd[ncol*itor+3]): continue
            iatoms = [ abs(amb_atoms_hyd[ncol*itor+i])/3 + 1
                    for i in range(ncol-1) ]
            index = amb_atoms_hyd[ncol*(itor+1) - 1]

            four_atoms.append( iatoms )
            force_consts.append( amb_consts[index-1])
            num_torsions.append( 1 )
            num_freqs.append( amb_freqs[index-1] )
            initial_phases.append( amb_phases[index-1] )

        # for torsion without hydrogen
        for itor in range( len(amb_atoms_woh)/ncol ):
            if not if_fun(amb_atoms_woh[ncol*itor+3]): continue
            iatoms = [ abs(amb_atoms_woh[ncol*itor+i])/3 + 1
                    for i in range(ncol-1) ]
            index = amb_atoms_woh[ncol*(itor+1) - 1]

            four_atoms.append( iatoms )
            force_consts.append( amb_consts[index-1])
            num_torsions.append( 1 )
            num_freqs.append( amb_freqs[index-1] )
            initial_phases.append( amb_phases[index-1] )

        all_torsion_info = dict(
                four_atoms=four_atoms,
                num_torsions=num_torsions,
                num_freqs=num_freqs,
                force_consts=force_consts,
                initial_phases=initial_phases )

        return all_torsion_info

    def _convert_coulomb(self):

        tpl = self.get_topology()
        # divede by 18.2223 to convert to charge in units of the electron charge
        amb_charges = [c/18.2223 for c in tpl.get_param_infos('CHARGE')]
        self.__coulomb_info = dict(charges=amb_charges)

    def _convert_vdw(self):

        tpl = self.get_topology()
        amb_types   = tpl.get_param_infos('ATOM_TYPE_INDEX')
        amb_indexes = tpl.get_param_infos('NONBONDED_PARM_INDEX')
        amb_c12s    = tpl.get_param_infos('LENNARD_JONES_ACOEF')
        amb_c6s     = tpl.get_param_infos('LENNARD_JONES_BCOEF')
        import math
        ntype = int(math.sqrt(len(amb_indexes)))

        c6s  = numpy.zeros((ntype, ntype))
        c12s = numpy.zeros((ntype, ntype))

        for itype_1 in range(ntype):
            itype = itype_1 + 1
            for jtype_1 in range(ntype):
                jtype = jtype_1 + 1

                ij = ntype*(itype-1) + jtype
                index = amb_indexes[ij-1]
                # if index < 0: print('negative', index, amb_c6s[index-1])
                if index < 0: continue
                c6s[itype_1, jtype_1]  = amb_c6s[index-1]
                c12s[itype_1, jtype_1] = amb_c12s[index-1]

        self.__vdw_info = dict(atom_types=amb_types, c6s=c6s, c12s=c12s)

if __name__ == '__main__':
    logger.log_level = logger.INFO
    # tpl = TopologyParser('./test/ala3.prmtop')
    fp = "./test/b2ar-act.prmtop.gz"
    
    tpl = TopologyParser(fp)
    tpl.parse()

    converter = Format2AmberBaseConverter(tpl)
    converter.convert()

    # natom_target = 1467
    # target_atoms = range(1, natom_target+1)
    # converter.apply_target_atoms(target_atoms)

    converter.print_atom()
    # converter.print_residue()
    # converter.print_bond()
    # converter.print_angle()
    # converter.print_torsion()
    # converter.print_improper()
    # converter.print_coulomb()
    # converter.print_vdw()
    # converter.print_bonded14_pairs()
    # converter.print_bonded_pairs()

    # converter.print_ibnd_to_ipair()
    # converter.print_iang_to_ipair()
    # converter.print_itor_to_ipair()
    # converter.print_iimp_to_ipair()
    # converter.print_i14_to_ipair()

    # converter.print_nonbonded_table()
