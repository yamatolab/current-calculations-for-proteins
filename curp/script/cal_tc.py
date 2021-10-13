"""Calculate autocorrelation function using flux data
CURP 1.1: Ishikura, 2016. Most of the Layout
CURP 1.2: Yamato, 2019. Heat Flux gestion (and commenting code).

- The type of flux data, i.e., scalar or vector, is specified by
  the newly added argument "-no_axes". With (without) this option,
  scalar (vector) flux is handled, and fortran module lib_acf (lib_hfacf)
  is employed. The lib_hfacf fortran module was newly added.
- For heat-flow, the unit of acf was should be changed to (A*(kcal/mol)/fs)^2.
  (to be worked)
"""

from __future__ import print_function

import sys
import time

import numpy

import netCDF4 as netcdf

# import parallel module
from curp import ParallelProcessor

# fortran module
from curp.script.lib_acf import cal_acf
from curp.script.lib_hfacf import cal_hfacf


def get_stringnames(string_array):
    return [''.join(string.tolist()).strip() for string in string_array]


def get_dt(flux_fn):
    """Get dt between each frame"""

    ncfile = netcdf.Dataset(flux_fn, mode='r')
    times = ncfile.variables['time'][:]
    d_t = times[1] - times[0]
    ncfile.close()

    return d_t


def gen_fluxdata(flux_fn, no_axes):

    t_0 = time.time()

    ncfile = netcdf.Dataset(flux_fn, mode='r')

    icom = 0
    donors = get_stringnames(ncfile.variables['donors'][:])
    acceptors = get_stringnames(ncfile.variables['acceptors'][:])
    npair = len(ncfile.dimensions['npair'])
    t_1 = time.time()
    print('preparing time: {:.3f} [s]'.format(t_1-t_0))

    ncfile.close()

    tsum = 0.0
    for ipair_1, (don, acc) in enumerate(zip(donors, acceptors)):
        print('loading for {}/{} {} {} ... ** '
              .format(ipair_1+1, npair, don, acc), end='')
        sys.stdout.flush()

        t_1 = time.time()
        ncfile = netcdf.Dataset(flux_fn, mode='r')

        if no_axes:
            flux = ncfile.variables['flux'][:, ipair_1, icom]
        else:
            flux = ncfile.variables['flux'][:, :, ipair_1, icom]

        ncfile.close()

        t_2 = time.time()
        print('load time: {:.3f} [s]'.format(t_2-t_1))
        tsum += t_2-t_1

        yield don, acc, flux

    print('total load time: {:.3f} [s]'.format(tsum))


class TCCalculator:
    def __init__(self, frame_range, avg_shift,
                 nsample, coef, no_axes, d_t=0.01):

        self.__fst_lst_intvl = frame_range
        fst, lst, intvl = self.__fst_lst_intvl
        self.nframe_acf = (lst - fst)/intvl + 1
        self.__fst = fst
        self.__lst = lst
        self.__intvl = intvl
        self.__shift = avg_shift
        self.__nsample = nsample
        self.__coef = coef
        self.d_t = d_t
        self.no_axes = no_axes

    def run(self, don, acc, flux):

        if self.no_axes:
            acf = cal_acf(flux, self.nframe_acf,
                                  self.__fst, self.__lst, self.__intvl,
                                  self.__shift, False, self.__nsample)
        else:
            acf = cal_hfacf(flux, self.nframe_acf,
                                      self.__fst, self.__lst, self.__intvl,
                                      self.__shift, False, self.__nsample, 3)

        tc = self.cal_tc(acf)
        return tc, acf

    def run_mpi(self, data, **other):
        t_0 = time.time()
        don, acc, flux = data

        if self.no_axes:
            acf = cal_acf(flux, self.nframe_acf,
                                  self.__fst, self.__lst, self.__intvl,
                                  self.__shift, False, self.__nsample)
        else:
            acf = cal_hfacf(flux, self.nframe_acf,
                                      self.__fst, self.__lst, self.__intvl,
                                      self.__shift, False, self.__nsample, 3)

        tc = self.cal_tc(acf)

        # self.print('cal time:', time.time()-t_0)
        print('    cal time: {:.3f} [s] for {} {}'
              .format(time.time()-t_0, don, acc))
        return don, acc, tc, acf

    def run_all(self, flux_iter):
        fst, lst, intvl = self.__fst_lst_intvl
        for flux in flux_iter:
            t_0 = time.time()

            if self.no_axes:
                acf = cal_acf(flux, self.nframe_acf, fst, lst, intvl,
                                      self.__shift, False, self.__nsample)
            else:
                acf = cal_hfacf(flux, self.nframe_acf, fst, lst,
                                          intvl, self.__shift, False,
                                          self.__nsample, 3)

            tc = cal_tc(acf, self.d_t, self.__coef)

            t_1 = time.time()
            print('    cal time: {:.3f} [s]'.format(t_1-t_0))

            yield tc, acf

    def cal_tc(self, acf):
        """
        Calculate final transport  coefficient.
        old version:  return self.__coef * self.d_t*1000.0 * numpy.sum(acf)
        """
        return self.__coef * self.d_t*1000.0 * numpy.trapz(acf,axis=0)[0]

    def get_times(self):
        return numpy.arange(0.0, self.nframe_acf*self.d_t, self.d_t)

    def print(self, *args, **kwds):
        print(*args, **kwds)


class TCWriter:

    """
    Write transport coefficient.
    """

    def __init__(self, fn):
        self.__fn = fn
        self.__file = open(self.__fn, 'wb')

    def write(self, don, acc, tc):
        fd = self.open()
        fd.write('{:>12} {:>12}  {}\n'.format(don, acc, tc))
        fd.flush()

    def open(self):
        if self.__file is None:
            fd = open(self.__fn, 'rb+')
            self.__file = fd

        return self.__file

    def close(self):
        if self.__file:
            self.__file.close()
            self.__file = None


class WriterBase:
    version = 0.1
    name = ''
    unit = ''

    def __init__(self, fp, nframe, d_t, label='test', title=''):

        self.__fp = fp
        self.__ncfile = None
        self.__nframe = nframe
        self.__title = title
        self.__label = label
        self.__d_t = d_t

        self.setup()

    def setup(self):

        ncfile = netcdf.Dataset(self.__fp, clobber=True,
                                mode='w', format='NETCDF3_64BIT')

        # set global attributes
        ncfile.title = self.__title
        ncfile.application = 'the CURP program'
        ncfile.program = 'cal-tc'
        ncfile.programVersion = str(0.7)
        ncfile.Convetsions = 'CURP'
        ncfile.ConvetsionVersion = str(self.version)

        # create dimensions
        ncfile.createDimension('npair', None)
        ncfile.createDimension('nframe', self.__nframe)
        ncfile.createDimension('nchar', 20)

        # create variables
        # time
        nc_time = ncfile.createVariable('time', 'f4', ('nframe',))
        nc_time.units = 'picosecond'
        d_t = self.__d_t
        times = [d_t * ifrm_1 for ifrm_1 in range(self.__nframe)]
        nc_time[:] = times[:]

        # donor and acceptor
        ncfile.createVariable('donors', 'c', ('npair', 'nchar'))
        ncfile.createVariable('acceptors', 'c', ('npair', 'nchar'))

        # trajectory
        nc_data = ncfile.createVariable(self.name, 'f4', ('npair', 'nframe'),)
        nc_data.units = self.unit

        ncfile.sync()
        ncfile.close()

    def write(self, ipair_1, don, acc, data):

        ncfile = self.open()

        nc_don = ncfile.variables['donors']
        nc_acc = ncfile.variables['acceptors']
        nc_data = ncfile.variables[self.name]

        nc_don[ipair_1] = list(don.ljust(20))
        nc_acc[ipair_1] = list(acc.ljust(20))
        nc_data[ipair_1] = data.ravel()

        self.close()

    def write_all(self, donors, acceptors, datas):

        ncfile = self.open()

        nc_don = ncfile.variables['donors']
        nc_acc = ncfile.variables['acceptors']
        nc_data = ncfile.variables[self.name]

        for ipair_1, (don, acc) in enumerate(zip(donors, acceptors)):
            nc_don[ipair_1] = list(don.ljust(20))
            nc_acc[ipair_1] = list(acc.ljust(20))

        nc_data[:] = datas

        self.close()

    def open(self):
        if self.__ncfile is None:
            ncfile = netcdf.Dataset(self.__fp, mode='r+')
            self.__ncfile = ncfile

        return self.__ncfile

    def close(self):
        if self.__ncfile:
            self.__ncfile.sync()
            self.__ncfile.close()
            self.__ncfile = None


class ACFWriter(WriterBase):
    name = 'acf'
    unit = '(kcal/mol/fs)^2'


def cal_tc(flux_fn, tc_fn="", acf_fn="", acf_fmt="netcdf",
           frame_range=[1, -1, 1], avg_shift=1, nsample=0, d_t=None,
           coef=1.0, use_debug=False, no_axes=False, **kwds):

    ti = time.time()
    par = ParallelProcessor()

    # determine delta t

    if d_t is None:
        d_t = get_dt(flux_fn)
        d_t = d_t * frame_range[2]    # in ps
    else:
        d_t = d_t * frame_range[2]    # in ps

    # prepare calculator
    cal = TCCalculator(frame_range, avg_shift, nsample, coef, no_axes, d_t)
    cal.print = par.write

    # prepare writer
    if par.is_root():
        tc_writer = TCWriter(tc_fn)

    # create ACF writer if it is necessary.
    if acf_fn:
        if par.is_root():
            acf_writer = ACFWriter(acf_fn, cal.nframe_acf, d_t)
    else:
        acf_writer = None

    # prepare flux data
    if par.is_root():
        fluxdata_iter = gen_fluxdata(flux_fn, no_axes)
    else:
        fluxdata_iter = None

    results_iter = par.run(cal.run_mpi, data=fluxdata_iter)

    for ipair_1, (don, acc, tc, acf) in enumerate(results_iter):
        if par.is_root():
            tc_writer.write(don, acc, tc)

            if acf_writer:
                t_acf = time.time()
                acf_writer.write(ipair_1, don, acc, acf)
                print('    writing acf: {:.4f} [s]'.format(time.time()-t_acf))

    if par.is_root():
        tc_writer.close()

    if par.is_root():
        print('total time: {:.3f} [s]'.format(time.time()-ti))

        print('==========================')
        print('    finished completely   ')
        print('==========================')



if __name__ == '__main__':
    from console import arg_cal_tc, exec_command

    parser = arg_cal_tc()
    exec_command(parser)
