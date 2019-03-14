
import os, sys
import traceback
import itertools as it
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)

################################################################################
from exception import CurpException

class FormatNotFoundError(CurpException): pass
class WriterNotFoundError(CurpException): pass
class NotFoundDictionary(CurpException): pass

################################################################################
class Writer:

    def __init__(self, trj_fp, trj_fmt,
            dt=0.01, is_vel=False, fst_lst_int=(1,-1,1)):

        fst, lst, intvl = fst_lst_int
        lst = 10**18 if lst == -1 else lst
        self.__fst_lst_int = (fst, lst, intvl)

        self.dt = dt

        if is_vel:
            Writer = _get_writer(trj_fmt, 'velocity')
        else:
            Writer = _get_writer(trj_fmt, 'coordinate')

        self.__writer = Writer(trj_fp, )
        self.__cnt = 0

    def write_header(self):
        self.__writer.write_header()

    def write(self, istp_1, trj, box=None):

        istp = istp_1+1
        fst, lst, intvl = self.__fst_lst_int

        if istp < fst: return
        if istp > lst: return

        if istp == fst:
            self.__writer.write(istp, trj, istp*self.dt, box)
            self.__cnt = 0
            return

        # if fst < istp < lst
        self.__cnt += 1

        if self.__cnt == intvl:
            self.__writer.write(istp, trj, istp*self.dt, box)
            self.__cnt = 0

    def close(self):
        self.__writer.close()


def _get_writer(format, writer_name):
    dictname = writer_name + '_writer_dict'

    # get the parsar's directory
    writer_path = os.path.dirname(__file__)

    # get the plugin names
    modpaths = []
    for fd in os.listdir(writer_path):
        abspath = os.path.join(writer_path, fd)
        if os.path.isdir(abspath) and not fd.startswith('_'):
            modpaths.append(abspath)

    for mpath in modpaths:
        modname = os.path.split(mpath)[-1]
        try:
            module = load_module(modname, writer_path)
            if hasattr(module, dictname):
                dictionary = getattr(module, dictname)
            else:
                raise NotFoundDictionary(dictname)

            if format in dictionary:
                Writer = dictionary[format]
                break

        except ImportError:
            print(traceback.format_exc())
            # raise FormatNotFoundError(format)

        except:
            raise

    else:
        raise WriterNotFoundError(format)

    return Writer

# class TrajectoryWriter:

    # def __init__(self, base_fp, writer, rotate=1, digit=3)
        # self.__base_fp = base_fp
        # self.__rotate  = rotate
        # self.__writer  = writer
        # self.__digit   = digit

        # self.__irot = 0


    # def write(self, snap):
        # self.__body.write(snap)

import imp
def load_module(modname,  basepath):
    f,n,d = imp.find_module(modname, [basepath])
    return imp.load_module(modname, f, n, d)


if __name__ == '__main__':

    Writer = _get_coordinate_writer('netcdf')
    print(Writer)


