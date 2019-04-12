#! /usr/bin/env python

"""Calculate autocorrelation function using flux data
CURP 1.1: Ishikura, 2016. Most of the Layout
CURP 1.2: Yamato, 2019. Heat Flux gestion (and commenting code).

- This task needs neither main() of curp.py nor .cfg file.
- The type of flux data, i.e., scalar or vector, is specified by
  the newly added argument "-no_axes". With (without) this option,
  scalar (vector) flux is handled, and fortran module lib_acf (lib_hfacf)
  is employed. The lib_hfacf fortran module was newly added.
- For heat-flow, the unit of acf was should be changed to (A*(kcal/mol)/fs)^2.
  (to be worked)
- In this update from 1.1 to 1.2, the following files are modified.
  $CURP_HOME/script/Makefile
  $CURP_HOME/script/cal_tc.py (this file).
"""

from __future__ import print_function
import os, sys
import math
import numpy
import gzip
import time

# import parallel module
from curp_module import ParallelProcessor

# fortran module
import lib_acf
import lib_hfacf

def parse_options():
    """Parse and get the command line options."""

    # make argument parser
    import argparse
    parser = argparse.ArgumentParser(description=(
        'Calculate transport coefficients from energy flux datas.') )

    # add argument definitions
    parser.add_argument(dest='flux_fn', metavar='FLUX_FILENAME',
            help=('The filepath of flux data.'))

    parser.add_argument('-o', '--tc-file', metavar='FILE',
            dest='tc_fn', required=True,
            help=('The filename of tc data.'))

    parser.add_argument('-a', '--acf-file', metavar='FILE',
            dest='acf_fn', required=False,
            default='',
            help=('The filename of acf data.'))

    parser.add_argument('-af', '--acf-file-format', metavar='FORMAT',
            dest='acf_fmt', required=False,
            default="netcdf", choices=['ascii', 'netcdf'],
            help=('The file format of acf data.'))

    parser.add_argument(
            '-r','--frame-range',metavar=('FIRST','LAST','INTER'),
            dest='frame_range', required=False,
            default=[1,-1,1], type=int, nargs=3,
            help='The range of frame; first frame, last frame, interval step.')

    parser.add_argument('-s', '--average-shift', metavar='SHIFT',
            dest='avg_shift', required=False,
            default=1, type=int,
            help='The frame to shift for averaging.')

    parser.add_argument(
            '--sample-number', dest='nsample',
            type=int, required=False, default=0,
            help=('number of sample for one flux data file. '
                'Default value is 0 that present to make samples '
                'as much as possible'))

    parser.add_argument('-dt', '--dt', metavar='DT',
            dest='dt', required=False,
            default=None, type=float,
            help=('t of between the neighbour frames. The unit is in ps. '
                + 'Default value is determined by time variable in flux data.'))

    parser.add_argument('-c', '--coefficient', metavar='COEFFICIENT',
            dest='coef', required=False,
            default=1.0, type=float,
            help='Multiply acf by given coefficient.')

    parser.add_argument(
            '-v', '--vervose', dest='use_debug',
            action='store_true', required=False,
            help='turn on debug mode.')

    parser.add_argument('-no_axes', dest='no_axes',
            action='store_true', required=False,
            help='with this option, scalar flux is handled.')

    # make arguments
    return parser.parse_args()

def get_stringnames(string_array):
    return [ ''.join(string.tolist()).strip() for string in string_array ]

def get_dt(flux_fn):

    import netCDF4 as netcdf
    ncfile = netcdf.Dataset(flux_fn, mode='r')
    times = ncfile.variables['time'][:]
    dt = times[1] - times[0]
    ncfile.close()

    return dt

def gen_fluxdata(flux_fn,no_axes):

    t0 = time.time()

    import netCDF4 as netcdf
    ncfile = netcdf.Dataset(flux_fn, mode='r')

    icom = 0
    donors    = get_stringnames (ncfile.variables['donors'][:] )
    acceptors = get_stringnames( ncfile.variables['acceptors'][:] )
    times     = ncfile.variables['time'][:]
    npair = len(ncfile.dimensions['npair'])
    t1 = time.time()
    print('preparing time: {:.3f} [s]'.format(t1-t0))

    ncfile.close()

    tsum = 0.0
    for ipair_1, (don,acc) in enumerate(zip(donors,acceptors)):
        print('loading for {}/{} {} {} ... ** '
                .format(ipair_1+1, npair, don, acc), end='')
        sys.stdout.flush()

        t1 = time.time()
        ncfile = netcdf.Dataset(flux_fn, mode='r')

        if no_axes:
            flux = ncfile.variables['flux'][:, ipair_1, icom]
        else:
            flux = ncfile.variables['flux'][:, :, ipair_1, icom]
        
        ncfile.close()

        t2 = time.time()
        print('load time: {:.3f} [s]'.format(t2-t1))
        tsum += t2-t1

        yield don,acc,flux

    print('total load time: {:.3f} [s]'.format(tsum))


class TCCalculator:

    def __init__(self, opts, dt=0.01):
        self.__fst_lst_intvl = opts.frame_range
        fst, lst, intvl = self.__fst_lst_intvl
        self.nframe_acf = (lst - fst)/intvl + 1
        self.__fst     = fst
        self.__lst     = lst
        self.__intvl   = intvl
        self.__shift   = opts.avg_shift
        self.__nsample = opts.nsample
        self.__coef    = opts.coef
        self.dt        = dt
        self.no_axes   = opts.no_axes

    def run(self, don, acc, flux):

        if self.no_axes: 
          acf = lib_acf.cal_acf(flux, self.nframe_acf,
                  self.__fst, self.__lst, self.__intvl, self.__shift,
                  False, self.__nsample)
        else:
          acf = lib_hfacf.cal_hfacf(flux, self.nframe_acf,
                  self.__fst, self.__lst, self.__intvl, self.__shift,
                  False, self.__nsample, 3)

        tc  = self.cal_tc(acf)
        return tc, acf

    def run_mpi(self, data, **other):
        t0 = time.time()
        don, acc, flux = data

        if self.no_axes:
          acf = lib_acf.cal_acf(flux, self.nframe_acf,
                  self.__fst, self.__lst, self.__intvl, self.__shift,
                  False, self.__nsample)
        else:
          acf = lib_hfacf.cal_hfacf(flux, self.nframe_acf,
                  self.__fst, self.__lst, self.__intvl, self.__shift,
                  False, self.__nsample, 3)

        tc  = self.cal_tc(acf)

        # self.print('cal time:', time.time()-t0)
        print('    cal time: {:.3f} [s] for {} {}'
                .format(time.time()-t0, don, acc))
        return don, acc, tc, acf

    def run_all(self, flux_iter):
        fst, lst, intvl = self.__fst_lst_intvl
        for flux in flux_iter:
            t0 = time.time()

            if self.no_axes:
              acf = lib_acf.cal_acf(flux, self.nframe_acf, fst, lst, intvl,
                      self.__shift, False, self.__nsample)
            else:
              acf = lib_hfacf.cal_hfacf(flux, self.nframe_acf, fst, lst, intvl,
                      self.__shift, False, self.__nsample, 3)

            tc  = cal_tc(acf, self.dt, self.__coef)

            t1 = time.time()
            print('    cal time: {:.3f} [s]'.format(t1-t0))

            yield tc, acf

    def cal_tc(self, acf):
        """Calculate final transport  coefficient."""
        return self.__coef * self.dt*1000.0 * numpy.sum(acf)

    def get_times(self):
        return numpy.arange(0.0, self.nframe_acf*self.dt, self.dt)

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


class Log:
    def __init__(self, filename):
        self.__fn = filename

        # initialize log file
        with open(filename, 'wb') as file:
            pass

    def write(self, line):
        with open(self.__fn, 'ab') as file:
            file.write(str(line) + '\n')

import netCDF4 as netcdf
class WriterBase:

    """
    """

    version = 0.1
    name = ''
    unit = ''

    def __init__(self, fp, nframe, dt, label='test', title='', compress=False):

        self.__fp     = fp
        self.__ncfile = None
        self.__nframe = nframe
        self.__title  = title
        self.__label  = label
        self.__dt     = dt

        self.setup()

    def setup(self):

        ncfile = netcdf.Dataset(self.__fp, clobber=True,
                mode='w', format='NETCDF3_64BIT')

        # set global attributes
        ncfile.title             = self.__title
        ncfile.application       = 'the CURP program'
        ncfile.program           = 'cal-tc'
        ncfile.programVersion    = str(0.7)
        ncfile.Convetsions       = 'CURP'
        ncfile.ConvetsionVersion = str(self.version)

        # create dimensions
        ncfile.createDimension('npair',  None)
        ncfile.createDimension('nframe', self.__nframe)
        ncfile.createDimension('nchar', 20)

        # create variables
        # time
        nc_time = ncfile.createVariable('time', 'f4', ('nframe',))
        nc_time.units = 'picosecond'
        dt = self.__dt
        times =  [ dt * ifrm_1 for ifrm_1 in range(self.__nframe) ]
        nc_time[:] = times[:]

        # donor and acceptor
        nc_don = ncfile.createVariable('donors',    'c', ('npair', 'nchar'))
        nc_acc = ncfile.createVariable('acceptors', 'c', ('npair', 'nchar'))

        # trajectory
        nc_data = ncfile.createVariable(self.name, 'f4', ('npair', 'nframe'),)
                # zlib=self.use_zlib, complevel=self.complevel)
        nc_data.units = self.unit

        ncfile.sync()
        ncfile.close()

    def write(self, ipair_1, don, acc, data):

        ncfile = self.open()

        nc_don   = ncfile.variables['donors']
        nc_acc   = ncfile.variables['acceptors']
        nc_data  = ncfile.variables[self.name]

        nc_don[ipair_1] = list(don.ljust(20))
        nc_acc[ipair_1] = list(acc.ljust(20))
        nc_data[ipair_1] = data.ravel()

        self.close()

    def write_all(self, donors, acceptors, datas):

        ncfile = self.open()

        nc_don  = ncfile.variables['donors']
        nc_acc  = ncfile.variables['acceptors']
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


def main():

    import time
    ti = time.time()
    # parse command line options
    opts = parse_options()

    par = ParallelProcessor()

    # log = Log(opts.log_fn)
    if opts.use_debug:
        log.write(opts)

    # determine delta t
    dt = get_dt(opts.flux_fn)

    if opts.dt is None:
        dt = dt * opts.frame_range[2] # in ps
    else:
        dt = opts.dt * opts.frame_range[2] # in ps

    # prepare calculator
    cal = TCCalculator(opts, dt)
    cal.print = par.write

    # prepare writer
    if par.is_root():
        tc_writer = TCWriter(opts.tc_fn)

    t3 = time.time()

    # create ACF writer if it is necessary.
    if opts.acf_fn:
        if par.is_root(): acf_writer = ACFWriter(opts.acf_fn,cal.nframe_acf,dt)
    else:
        acf_writer = None

    # prepare flux data
    if par.is_root():
        fluxdata_iter = gen_fluxdata(opts.flux_fn,opts.no_axes)
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
    main()
