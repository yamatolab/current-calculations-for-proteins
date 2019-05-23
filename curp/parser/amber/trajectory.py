from __future__ import print_function

scale_factor = 20.455 * 10.0**(-3)

try:
    import lib_parser
except ImportError:
    lib_parser = None

from exception import CurpException
class NumAtomInvalidError(CurpException): pass


import gzip
import numpy

class ParserPyBase:

    def __init__(self, filename, natom, use_pbc=False):
        self.__natom = natom
        self.__use_pbc = use_pbc

        if filename.endswith('.gz'):
            self.__file = gzip.open(filename, 'rb')
        else:
            self.__file = open(filename, 'rb')
        title = self.parse_header(self.__file)

    def _set_module(self, module):
        pass

    def __next__(self):
        crdvel = list(self.parse_crdvel_iter(self.__file, self.__natom))
        box = self.parse_pbc(self.__file) if self.__use_pbc else None
        if len(crdvel) == 0:
            raise StopIteration()

        return numpy.array(crdvel), box

    next = __next__ # for python 2.x

    def __iter__(self):
        return self

    def parse_header(self, file):
        return file.next()

    def parse_data(self, file, natom):
        nline = natom*3 / 10
        rem   = natom*3 % 10
        try:
            for iline in range(nline):
                line = file.next()
                for icol in range(10):
                    beg, end = 8*icol, 8*(icol+1)
                    yield float(line[beg:end])

            else:
                if rem != 0:
                    line = file.next()
                    for icol in range(rem):
                        beg, end = 8*icol, 8*(icol+1)
                        yield float(line[beg:end])

        except ValueError:
            msg = "Given the number of atoms is {}.".format(natom)
            raise NumAtomInvalidError(msg)

    def parse_crdvel_iter(self, file, natom):
        gen_data = self.parse_data(file, natom)

        for iatm in range(natom):
            x = gen_data.next()
            y = gen_data.next()
            z = gen_data.next()
            yield x, y, z

    def parse_pbc(self, file):

        # parse the info of periodic boudary condition
        line = file.next()
        p1, p2, p3 = line[0:8],  line[8:16], line[16:24]
        return [float(p) for p in [p1, p2, p3]]

    def close(self):
        self.__file.close()

    def get_dt(self):
        return 0.1


class ParserFortBase:

    def __init__(self, filename, natom, use_pbc=False):
        self.__mod.initialize(filename)
        self.__natom = natom
        self.__use_pbc = use_pbc

    def _set_module(self, module):
        self.__mod = module

    def __next__(self):
        crd, ierr = self.__mod.parse_next(self.__natom)
        if ierr != 0:
            self.__mod.finalize()
            raise StopIteration()

        if self.__use_pbc:
            box, ierr = self.__mod.parse_pbc(self.__natom)
            return crd, box
        else:
            return crd, None
    next = __next__ # for python 2.x

    def __iter__(self):
        return self

ParserBase = ParserPyBase

class CoordinateParser(ParserBase):

    def __init__(self, *args, **kwds):
        if lib_parser:
            self._set_module(lib_parser.formatted_coordinate)
        ParserBase.__init__(self, *args, **kwds)

class VelocityParser_half_dt(ParserBase):

    def __init__(self, *args, **kwds):
        if lib_parser:
            self._set_module(lib_parser.formatted_velocity)
        ParserBase.__init__(self, *args, **kwds)

        self.__vel_old, box = ParserBase.__next__(self)

    def __next__(self):
        vel, box = ParserBase.__next__(self)
        vel_half = 0.5 * (vel + self.__vel_old)
        self.__vel_old = vel
        return scale_factor*vel_half, box

    next = __next__ # for python 2.x

class VelocityParserNormal(ParserBase):

    def __init__(self, *args, **kwds):
        if lib_parser:
            self._set_module(lib_parser.formatted_velocity)
        kwds['use_pbc'] = False
        ParserBase.__init__(self, *args, **kwds)

    def __next__(self):
        vel, box = ParserBase.__next__(self)
        return scale_factor * vel, box

    next = __next__ # for python 2.x

VelocityParser = VelocityParserNormal

import netCDF4 as netcdf

class NetCDFReaderBase:

    trjname = None # A inheritanced class must have this variable.
    scale_factor   = 1.0

    def __init__(self, filename, natom=None, use_pbc=False):

        self.__fp = filename

        # read NetCDF file
        try:
            from netCDF4 import Dataset
        except ImportError:
            raise
        trjobj = Dataset(filename, mode='r')

        self._natom  = len(trjobj.dimensions['atom'])
        self._nframe = len(trjobj.dimensions['frame'])
        self.title   = getattr(trjobj, 'title', None)
        if use_pbc:
            if 'cell_lengths' in trjobj.variables: #add 17.06.29
                self._box1 = trjobj.variables['cell_lengths']
            if 'cell_angles' in trjobj.variables:  #add 17.06.29
                self._box2 = trjobj.variables['cell_angles']
        self.__use_pbc = use_pbc

        self.__iter_cnt = 0

        trjobj.sync()
        trjobj.close()
        self.__trjobj = None

    def __getitem__(self, ifrm):
        trjobj = self.open()
        trj_vars = trjobj.variables
        time = trj_vars['time'][ifrm]
        crdvel = trj_vars[self.trjname][ifrm]

        if self.__use_pbc:
            if 'cell_lengths' in trj_vars: #add 17.06.29
                box1 = trj_vars['cell_lengths'][ifrm]
            else:  #add 17.06.29
                box1=''  #add 17.06.29
            if 'cell_angles' in trj_vars:  #add 17.06.29
                box2 = trj_vars['cell_angles'][ifrm]
            else:  #add 17.06.29
                box2=''  #add 17.06.29
        else:
            n = len(crdvel)
            box1, box2 = None, None

        self.close()

        #TODO convert the box information 1 matrix.
        box = box1

        return crdvel*self.scale_factor, box

    def __next__(self):
        try:
            crdvel, box = self[self.__iter_cnt]
            self.__iter_cnt += 1

        except IndexError:
            self.__iter_cnt = 0
            raise StopIteration()

        return crdvel, box

    next = __next__ # for python 2.x
    def __iter__(self): return self

    def parse_header(self, file):
        return file.next()

    def open(self):
        if self.__trjobj is None:
            trjobj = netcdf.Dataset(self.__fp, mode='r')
            self.__trjobj = trjobj

        return self.__trjobj

    def close(self):
        if self.__trjobj:
            self.__trjobj.sync()
            self.__trjobj.close()
            self.__trjobj = None

    def get_dt(self):
        trjobj = self.open()
        trj_vars = trjobj.variables
        times = trj_vars['time'][ifrm]
        dt = times[1] - times[0]
        self.close()
        return dt

class NetCDFReaderBaseNew:

    trjname = None # A inheritanced class must have this variable.
    scale_factor   = 1.0

    def __init__(self, filename, natom=None, use_pbc=False):

        self.__filename = filename
        trjobj = self.open(filename)

        self._natom  = len(trjobj.dimensions['atom'])
        self._nframe = len(trjobj.dimensions['frame'])
        self.title   = getattr(trjobj, 'title', None)
        if use_pbc:
            self._box1 = trjobj.variables['cell_lengths']
            self._box2 = trjobj.variables['cell_angles']
        self.__use_pbc = use_pbc

        self.__iter_cnt = 0

    def open(self, filename):
        # read NetCDF file
        try:
            from netCDF4 import Dataset
        except ImportError:
            raise
        trjobj = Dataset(filename)
        return trjobj


    def __getitem__(self, ifrm):
        trjobj = self.open(self.__filename)
        trj_vars = trjobj.variables
        time = trj_vars['time'][ifrm]
        crdvel = trj_vars[self.trjname][ifrm]

        if self.__use_pbc:
            box1 = trj_vars['cell_lengths'][ifrm]
            box2 = trj_vars['cell_angles'][ifrm]
        else:
            n = len(crdvel)
            box1, box2 = None, None

        trjobj.close()

        return crdvel*self.scale_factor, (box1, box2)

    def __next__(self):
        try:
            crdvel, box = self[self.__iter_cnt]
            self.__iter_cnt += 1

        except IndexError:
            self.__iter_cnt = 0
            raise StopIteration()

        return crdvel, box

    next = __next__ # for python 2.x
    def __iter__(self): return self

    def parse_header(self, file):
        return file.next()

    def close(self):
        pass
        # self.__trjobj.close()


class NetCDFCoordinateReader(NetCDFReaderBase):

    trjname = 'coordinates'

    def __init__(self, *args, **kwds):
        NetCDFReaderBase.__init__(self, *args, **kwds)


class NetCDFVelocityReader(NetCDFReaderBase):

    trjname = 'velocities'
    # The velocity unit in NetCDF is angstrom/(ps/20.455)
    scale_factor = 1.0/10.0**3  # A/ps <= A/fs

    def __init__(self, *args, **kwds):
        NetCDFReaderBase.__init__(self, *args, **kwds)


class InvalidAtomNumber: pass
class RestartParser:

    def __init__(self, filename, natom, use_pbc=False, trj_type='crdvel'):
        self.__filename = filename
        self.__natom    = natom
        self.__use_pbc  = use_pbc

    @classmethod
    def set_trjtype(cls, trj_type='crd'):

        cls.__use_crd = False
        cls.__use_vel = False
        cls.__parsed = False

        if trj_type == 'crd':
            cls.__use_crd = True
        elif trj_type == 'vel':
            cls.__use_vel = True
        else:
            pass

    def __next__(self):

        if self.__parsed:
            raise StopIteration()

        crd, vel, pbc = self.parse()

        if self.__use_crd:
            self.__parsed = True
            return crd, pbc
        elif self.__use_vel:
            self.__parsed = True
            return vel, pbc
        else: 
            self.__parsed = True
            return (crd, vel), pbc

    next = __next__ # for python 2.x

    def __iter__(self):
        return self

    def parse(self):
        with open(self.__filename, 'rb') as file:
            title, natom, simtime = self.parse_header(file)
            import numpy
            crd = list( self.parse_crd_iter(file, natom) )
            vel = list( self.parse_vel_iter(file, natom) )

            if self.__use_pbc:
                pbc = self.parse_pbc(file)
            else:
                pbc = None

        return numpy.array(crd), numpy.array(vel), pbc

    def close(self):
        pass

    def parse_header(self, file):
        # parse title line
        title = file.next()

        # parse info line
        line = file.next()
        cols = line.split()
        natom, simtime = int(cols[0]), float(cols[1])

        if natom != self.__natom:
            raise InvalidAtomNumber()

        return title, natom, simtime

    def parse_data(self, file, natom):
        for iline in range(natom/2):
            line = file.next()
            for icol in range(6):
                beg, end = 12*icol, 12*(icol+1)
                yield float(line[beg:end])

        else:
            if natom%2 == 1:
                line = file.next()
                for icol in range(3):
                    beg, end = 12*icol, 12*(icol+1)
                    yield float(line[beg:end])

    def parse_crd_iter(self, file, natom):
        gen_data = self.parse_data(file, natom)

        for iatm in range(natom):
            x = gen_data.next()
            y = gen_data.next()
            z = gen_data.next()
            yield x, y, z

    def parse_vel_iter(self, file, natom):
        gen_data = self.parse_data(file, natom)
        c = scale_factor

        for iatm in range(natom):
            x = gen_data.next() * c
            y = gen_data.next() * c
            z = gen_data.next() * c
            yield x, y, z

    def parse_pbc(self, file):

        # parse the info of periodic boudary condition
        try:
            line = file.next()
            p1, p2, p3 = line[0:12],  line[12:24], line[24:36]
            p4, p5, p6 = line[36:48], line[48:60], line[60:72]
            pbc = [float(p) for p in [p1, p2, p3, p4, p5, p6]]
        except StopIteration:
            pbc = None

        return None


class TrajectoryWriterBase:

    suffix = ''

    def __init__(self, filepath, fst_lst_int=(1,-1,1),
            use_pbc=False,rotate=1,compress=0):
        self.__fp = filepath
        self.__rotate = rotate
        self.__fst_lst_int = fst_lst_int
        self.__compress = compress
        self.__use_pbc = use_pbc

    def write_header(self):
        pass

    def write_all(self, trj):
        for snap in trj:
            self.write(snap)

        self.close()

    def open(self):
        fp = self.get_fp()
        self.__file = open(fp, 'wb')
        return self.__file

    def write(self, ifrm, snap, time, box):
        pass

    def check(self):
        fp = self.get_fp()
        if os.path.exists(fp):
            raise TrajectoryExistsError(fp)

    def close(self):
        self.__file.close()

    def get_fp(self):
        return self.__fp

    def use_pbc(self):
        return self.__use_pbc

import lib_trj
class AsciiCoordinateWriter(TrajectoryWriterBase):

    suffix = 'mdcrd'

    def __init__(self, *args, **kwds):
        TrajectoryWriterBase.__init__(self, *args, **kwds)
        self.write_header()

    def write_header(self):
        msg = "This trajectory was generated by the CURP."
        file = self.open()
        file.write(msg+'\n')
        file.close()

    def write(self, ifrm, crd, time=None, box=None):
        fp = self.get_fp()

        # write next coordinate
        lib_trj.write_trajectory(crd, fp)


class AsciiVelocityWriter(TrajectoryWriterBase):

    suffix = 'mdvel'
    scale_factor = 10.0**3 / 20.455

    def __init__(self, *args, **kwds):
        TrajectoryWriterBase.__init__(self, *args, **kwds)
        self.write_header()

    def write_header(self):
        msg = "This trajectory was generated by the CURP."
        file = self.open()
        file.write(msg+'\n')
        file.close()

    def write(self, ifrm, vel, time=None, box=None):
        fp = self.get_fp()

        # write next velocity
        # convert curp unit into amber unit for velocity
        c = self.scale_factor
        lib_trj.write_trajectory(vel*c, fp)


class NetCDFWriterBase(TrajectoryWriterBase):

    # Please set these class variables in a inheritant class.
    suffix   = 'xxx.nc'
    trj_type = ''
    units    = ''
    ncdf_scale = 1.0
    scale_factor = 1.0

    # class variable
    version  = 1.0

    def __init__(self, *args, **kwds):
        self.__ncfile = None
        self.__ifrm   = 0
        TrajectoryWriterBase.__init__(self, *args, **kwds)

    def setup(self, natom):
        """Setup for NetCDF dataset."""

        ncfile = netcdf.Dataset(self.get_fp(), clobber=True,
                mode='w', format='NETCDF3_64BIT')

        # set global attributes
        ncfile.title             = 'netCDF '+self.trj_type +'trajectory'
        ncfile.application       = 'CURP'
        ncfile.program           = 'Utility tools in the CURP'
        ncfile.programVersion    = '0.6'
        ncfile.Convetsions       = 'AMBER'
        ncfile.ConvetsionVersion = str(self.version)

        # create dimensions
        ncfile.createDimension('frame', None)
        ncfile.createDimension('spatial', 3)
        ncfile.createDimension('atom', natom)
        ncfile.createDimension('cell_spatial', 3)
        ncfile.createDimension('cell_angular', 3)
        # file.createDimensions('label', ?)

        # create variables
        time = ncfile.createVariable('time', 'f4', ('frame',))
        time.units = 'picosecond'

        trj = ncfile.createVariable(self.trj_type, # coordinates or velocities
                'f4', ('frame', 'atom', 'spatial' ))
        trj.units = self.units
        trj.scale_factor = self.ncdf_scale


        cell_lengths = ncfile.createVariable('cell_lengths',
                'f8', ('frame', 'cell_spatial'))
        cell_lengths.units = 'angstrom'
        
        cell_angles = ncfile.createVariable('cell_angles',
                'f8', ('frame', 'cell_angular'))
        cell_angles.units = 'degree'

        # file.createVariable('spatial', 'c')
        # file.variables['spatial'] = 'xyz'
        # print(file.variables)

        # sync
        ncfile.sync()
        ncfile.close()

    def write_header(self):
        """Write header."""
        pass

    def write(self, ifrm, snap, time=None, box=None):
        """Write snapshots to NetCDF file object."""

        # setup at first step
        if self.__ifrm == 0:
            natom = len(snap)
            self.setup(natom)

        ncfile = self.open()
        # ncfile = self.__ncfile

        self.__ifrm += 1
        ifrm_1 = self.__ifrm - 1

        nc_trj = ncfile.variables[self.trj_type]
        nc_time = ncfile.variables['time']

        nc_trj[ifrm_1] = snap[:]/self.scale_factor
        nc_time[ifrm_1] = time

        if isinstance(box,numpy.ndarray): # periodic boundary conditios ##add 17.06.29
            print('here', box)
            ncfile.variables['cell_lengths'][ifrm_1] = box
            ncfile.variables['cell_angles'][ifrm_1]  = box
            # make this code right in next version.

        self.close()

    def open(self):
        if self.__ncfile is None:
            ncfile = netcdf.Dataset(self.get_fp(), clobber=True, mode='r+')
            self.__ncfile = ncfile

        return self.__ncfile

    def close(self):
        if self.__ncfile:
            self.__ncfile.sync()
            self.__ncfile.close()
            self.__ncfile = None
        

class NetCDFCoordinateWriter(NetCDFWriterBase):

    suffix   = 'crd.nc'
    trj_type = 'coordinates'
    units    = 'angstrom'

    def __init__(self, *args, **kwds):
        NetCDFWriterBase.__init__(self, *args, **kwds)


class NetCDFVelocityWriter(NetCDFWriterBase):

    suffix   = 'vel.nc'
    trj_type = 'velocities'
    units    = 'angstrom/picosecond'
    ncdf_scale = 20.455
    scale_factor = 1.0/10.0**3 # fs => ps

    def __init__(self, *args, **kwds):
        NetCDFWriterBase.__init__(self, *args, **kwds)


if __name__ == "__main__":

    import os

    curp_path = os.environ['CURP_HOME']
    exapmle_path = prmtop_fn = os.path.join(
            curp_path,'test','amber-enk-vacuum')

    prmtop_fn = os.path.join(exapmle_path, 'prmtop_fn')
    crds_fn = os.path.join(exapmle_path, 'sam.nccrd')
    vels_fn = os.path.join(exapmle_path, 'sam.ncvel')

    from benchmarker import Benchmarker
    with Benchmarker(width=20) as bm:
        # with bm('parse'):
            # rst_fn = "test/system.rst"
            # parser = RestartParser(rst_fn, 30661)
            # crd, vel, pbc = parser.parse()
            # crd, vel = parser.parse()
            # print(vel)
        # with bm('print'):
        #     crd, vel
        #     print()
        #     for p in crd:
        #         print(p)
        #     for p in vel:
        #         print(p)

        # print('crd')
        # for c in crd:
        #     print(c)
        # print('vel')
        # for v in vel:
        #     print(v)
        # print('pbc')
        # print(pbc)

        # with bm('crd'):
        #     print()
        #     parser = CoordinateParser('test/sam-nwat.mdcrd', 1799)
        #     i = 0
        #     for crd, box in parser:
        #         i += 1
        #         if i%100==0:
        #             print('{:>5} x 10^2: size = {:>5}'.format(
        #                 i/100, len(crd)))
        #     print(i)
                    
        # with bm('crd+pbc'):
            # print()
            # parser = CoordinateParser('test/ala3-woct.mdcrd.gz',
                    # 5844, use_pbc=True)
            # i = 0
            # for crd, box in parser:
                # i += 1
                # if i%100==0:
                    # print('{:>5} x 10^2: size = {:>5}'.format(
                        # i/100, len(crd)))

            # print(i)

        with bm('NetCDF'):
            crd_trj = NetCDFCoordinateReader(crds_fn)
            for crd, pbc in crd_trj:
                crd
                # print(crd)

            vel_trj = NetCDFVelocityReader(vels_fn)
            for vel, pbc in vel_trj:
                vel
                # print(vel/scale_factor) # A/fs => A/(ps/20.455)

