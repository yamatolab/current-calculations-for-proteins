from __future__ import print_function

import os, sys
from collections import OrderedDict as odict
import numpy

# curp mudules
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)
import exception
import clog as logger
import interact_table as it

class InvalidGroupName(exception.CurpException): pass

import ini_parser as ini
class GroupPairParser(ini.IniParser):

    def __init__(self, gpair_fn, gnames):
        ini.IniParser.__init__(self, gpair_fn)
        self.__gnames_base = gnames
        self.parse()
        self.check()

    def check(self):

        gnames = self.get_secnames()

        if not (set(gnames) <= set(self.__gnames_base)):
            invalid_names = list(set(gnames) - set(self.__gnames_base))
            msg = "{} in {} is invalid group names !".format(
                    invalid_names, self.get_filename())
            
            raise InvalidGroupName(msg)

        for gname_i in gnames:
            gnames_j = self[gname_i]

            if not (set(gnames_j) <= set(self.__gnames_base)):
                invalid_names = list(set(gnames_j) - set(self.__gnames_base))
                msg = "{} in {} is invalid group names !".format(
                        invalid_names, self.get_filename())
                raise InvalidGroupName(msg)

    def get_gpair_table(self):
        return list(ini.IniParser.gen_secname_cols_pair(self))


class GroupPair:

    def __init__(self, gpair_table, gname_to_iatoms, natom):
        self.__gpair_table = gpair_table
        self.__gname_to_iatoms = gname_to_iatoms
        self.__natom = natom
        self.__table_with_gpair = None

    def gen_inttable(self):
        """Generate the object equivalent with interaction table
        from group pair table, ex)
        """
        gname_to_iatoms = dict(self.__gname_to_iatoms)

        for gname_i, jgnames in self.__gpair_table:

            jatoms = []
            for gname_j in jgnames:
                jatoms += gname_to_iatoms[gname_j]

            comp_jatoms = list(gen_compressed_iatms(jatoms))
            for iatm in gname_to_iatoms[gname_i]:
                for jatm_beg, jatm_end in comp_jatoms:
                    yield iatm, jatm_beg, jatm_end

    def get_inttable_with_gpair(self, base_table):
        if self.__table_with_gpair is None:
            self.__table_with_gpair = self._make_table_with_gpair(base_table)
        return self.__table_with_gpair

    def _make_table_with_gpair(self, base_table):
        """Make interaction table with group pairs."""

        table_with_gpair = numpy.array(list(self.gen_inttable()))
        base_table = numpy.array(list(base_table))

        import lib_group_pair
        lib_gpair = lib_group_pair.within_gpair
        lib_gpair.setup(table_with_gpair, self.__natom)

        # make interaction table
        new_table, ntable = lib_gpair.get_nonbonded_table(base_table)
        new_table = new_table[:ntable].tolist()
        return it.InteractionTable(base_table=new_table)

def gen_compressed_iatms( iatoms ):
    """Convert: [1,2,3,6,7,8,...,100] => [(1,3),(6,100)]"""

    iatoms.sort()
    iatm_beg = iatoms[0]
    iatm_end = iatoms[0]

    for iatm in iatoms[1:]:

        if iatm != iatm_end+1:
            yield (iatm_beg, iatm_end)
            iatm_beg = iatm

        iatm_end = iatm

    else:
        yield (iatm_beg, iatm_end)


if __name__ == '__main__':
    import os
    gpair_fn = os.path.join(os.environ['CURP_HOME'],
            'test', 'energy-flux', 'pickup', 'gpair.dat')
    tplprm_fn = os.path.join(os.environ['CURP_HOME'],
            'test', 'amber-b2AR', 'stripped.prmtop.gz')

    from parser.amber.topology import TopologyParser, Format2AmberBaseConverter
    raw_tpl = TopologyParser(tplprm_fn)
    raw_tpl.parse()
    tpl = Format2AmberBaseConverter(raw_tpl, use_atomtype=False)
    tpl.convert()

    atom_info = tpl.get_atom_info()
    res_info = tpl.get_residue_info()

    import group
    gnames = [ gname for gname, iatoms
            in group.gen_residue_group(res_info, atom_info) ]

    # for gname in gnames:
        # print(gname)

    gpair_parser = GroupPairParser(gpair_fn, gnames)
    for i in gpair_parser.get_gpair_table():
        print(i)

    # with bm('gen_compressed_iatms'):
        # print()
        # iatoms = [1,2,3,8,9,10,11,12,13]
        # table = gen_compressed_iatms(iatoms)
        # print(list(table))

