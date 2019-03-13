from __future__ import print_function

import os, sys
from collections import OrderedDict as odict

topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)
import exception

class IniParser:

    """
    $ cat target_file
    [SEC_1]
    attr1 = val1
    attr2 = val2
    1 2 3 4
    [SEC_2]
    5 6 7 8 9 10
    .
    .
    .
    [SEC_N]
    attr1 = val1
    attr2 = val2
    100 101 102

    =>

    {'SEC1':[1,2,3,4], 'SEC_2' [5,6,7,8,9,10], ..., 'SEC_N': [100,101,102]}
    """

    section_end = "#[END_SECTION]"
    comment = '#'

    def __init__(self, filename):
        self.__filename = filename
        self.parse()

    def __getitem__(self, secname):
        return self.get_vals(secname)

    def get_filename(self):
        return self.__filename
        
    # def __setitem__(self, secname, value):
    #     pass

    def get_secnames(self):
        """Get the section names in order."""
        return self.__secname_to_vals.keys()

    def get_vals(self, secname):
        """Get values for section name."""
        return self.__secname_to_vals[secname]

    def get_attr_dict(self, secname):
        """Get attribute value dictionary for section name."""
        return self.__secname_to_attrdict[secname]

    def parse(self):
        """Parse ini format content."""

        with open(self.__filename, 'r') as file:
            secname_lines_pairs = self._parse_secname_lines_pairs(file)

        # parse vals
        self.__secname_to_vals = odict()
        for secname, lines in secname_lines_pairs:
            self.__secname_to_vals[secname] = list(self._gen_col(lines))

        # parse attributes
        self.__secname_to_attrdict = odict()
        for secname, lines, in secname_lines_pairs:
            self.__secname_to_attrdict[secname] = odict(
                    self._gen_attr_val(lines))

    def gen_secname_cols_pair(self):
        """Generate section_name: columns pair."""
        for secname in self.get_secnames():
            yield secname, self[secname]

    def _gen_col(self, lines):
        for line in lines:
            if '=' in line: continue
            for col in line.split():
                yield col

    def _gen_attr_val(self, lines):
        for line in lines:
            cols = line.split('=')
            if len(cols) == 2:
                yield cols[0].strip(), cols[1].strip()

    def _parse_secname_lines_pairs(self, file):
        sections = []
        gen_lines = self._gen_optimized_lines(file)

        for line in gen_lines:
            if self._is_section(line):
                lines = list(self._gen_section_lines(gen_lines))
                sections.append( (line[1:-1].strip(), lines) )

        return sections

    def _gen_section_lines(self, gen):
        for line in gen:
            if line == self.section_end:
                break
            yield line

    def _gen_optimized_lines(self, file):
        """Generate lines with the end of sections and without the comments."""

        # search first section marker
        for line in file:
            line = line.split(self.comment)[0].strip()
            if self._is_section(line):
                yield line
                break
            elif line == '':
                continue
            else:
                mes = 'found the invalid string before first section marker: '
                raise ParserError(mes + line)

        # rest
        for line in file:
            line = line.split(self.comment)[0].strip()
            if self._is_section(line):
                yield self.section_end
                yield line
            elif line != '':
                yield line
            else:
                continue
        else:
            yield self.section_end

    def _is_section(self, line):
        if line == self.section_end:
            return False
        elif line.startswith('[') and line.endswith(']'):
            return True
        elif line.startswith('[') or line.endswith(']'):
            raise SectionNotFoundError(line)
        else:
            return False


