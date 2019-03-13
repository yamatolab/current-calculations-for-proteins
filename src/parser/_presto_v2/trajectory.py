import os, sys
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)

from exception import CurpException

class FileNotReadError(CurpException): pass

import lib_parser
class ParserBase:

    def __init__(self, filename, natom):
        self.__mod.initialize(filename)
        self.__natom = natom

    def _set_module(self, module):
        self.__mod = module

    def __next__(self):
        crd, ierr = self.__mod.parse_next(self.__natom)
        if ierr != 0:
            self.__mod.finalize()
            raise StopIteration()
        else:
            return crd
    next = __next__ # for python 2.x

    def __iter__(self):
        return self


class CoordinateParser(ParserBase):

    def __init__(self, filename, natom):
        self._set_module(lib_parser.unformatted_coordinate)
        ParserBase.__init__(self, filename, natom)


class VelocityParser(ParserBase):

    def __init__(self, filename, natom):
        self._set_module(lib_parser.unformatted_velocity)
        ParserBase.__init__(self, filename, natom)


class RestartParser:

    def __init__(self, filename, natom):
        self.__mod = lib_parser.unformatted_restart
        self.__mod.initialize(filename)
        self.__natom = natom
        self.__mod.parse1()

    def parse(self):
        crd, ierr = self.__mod.parse_crd(self.__natom)
        if ierr != 0:
            self.__mod.finalize()
            mes = 'in reading coodinate'
            raise FileNotReadError(mes)

        vel, ierr = self.__mod.parse_vel(self.__natom)
        if ierr != 0:
            self.__mod.finalize()
            mes = 'in reading velocity'
            raise FileNotReadError(mes)

        return crd, vel


if __name__ == "__main__":

    # import lib_parser

    # help(lib_parser)

    rst_fn = "./test/pyp_md.rst"

    parser = RestartParser(rst_fn, 10779)
    crd, vel = parser.parse()
    print( crd)
    print( vel)
