"""Write the results of curp calculator

CURP 1.1: Ishikura, 2016. Most of the Layout  
CURP 1.2: Laprevote, 2019. Heat Flux gestion (and commenting code).
CURP 1.2: Yamato, 2021. Minor modification for intra-residue heat flux.

In main() of curp.py, a writer is initialized along with the calculator.
The writer varies depending on setting options (in .cfg) contained in
a setting object.
The writer will vary because of (in this order):
 * Decided in get_writer:
  - The calculation method (momentum-current, energy or heat-flux)
 * Decided in MultiFluxWriter and MultiCurrentWriter:
  - The output file format (ascii or netcdf)
  - The use of groups or not
 * Decided in CurrentWriter, FluxWriter and NetCDFFluxWriter:
  - The compression of the output file
  - In FluxWriter and NetCDFFluxWriter, the presence of axes
    (as heat flux is a vector and energy flux a single value)

Each of the XxxWriter will format lines using the output of the
calculator for one frame, and use either Writer or AxisWriter object
to write them in a file.
"""

from __future__ import print_function

import os
import numpy


def get_writer(setting, decomp_list, target_anames, group_names,
               gpair_table=None):
    """Returns Writer object depending on the setting.curp.method."""

    method = setting.curp.method
    if method == 'momentum-current':
        obj = MultiCurrentWriter(
                setting=setting, decomp_list=decomp_list,
                target_anames=target_anames, group_names=group_names)

    elif method == 'energy-flux':
        obj = MultiFluxWriter(
                setting=setting, decomp_list=decomp_list,
                target_anames=target_anames, group_names=group_names,
                gpair_table=gpair_table)

    elif method == 'heat-flux':
        obj = MultiFluxWriter(
                setting=setting, decomp_list=decomp_list,
                target_anames=target_anames, group_names=group_names,
                gpair_table=gpair_table, axes=('i','j','k'))

    else:
        pass

    return obj


def write(setting, results, formats):
    key_to_writer = {}
    if setting.output.energy:
        key_to_writer['energy'] = Writer(
                energy.filename, formats['energy'],
                energy.frequency, energy.compress)

    if setting.output.force:
        key_to_writer['force'] = Writer(
                force.filename, formats['force'],
                force.frequency, force.compress)

    if setting.output.twobody_force:
        key_to_writer['twobody_force'] = Writer(
                twobody_force.filename, format['twobody_force'],
                twobody_force.frequency,twobody_forceforce.compress)

    for res in results:
        for key, writer in key_to_writer.items():
            data = res[key]
            writer.write(data)


class Writer:
    """Called to write files.

    Used in CurrentWriter.
    Parent of AxisWriter.

    Keyword arguments:
    filename -- string, will precede the name of written files
    frequency -- integer 
    compresslevel -- integer, ==0 when compress=no in .cfg file
    use_axes -- boolean, True when Writer is used in AxisWriter.
    """

    def __init__(self, filename, frequency,
                 compresslevel=6, use_axes=False):
        self._base_fn = filename
        self._nfreq = frequency
        self._ifreq = 0     # The number of the frame
        self._counter = 0

        # Choose the open function used
        self._use_zip = compresslevel != 0
        if self._use_zip:
            import gzip
            def gzopen(filename, mode):
                import gzip 
                return gzip.open(filename, mode, compresslevel)
            self._open = gzopen
        else:
            self._open = open

        self.make_filename(self._ifreq, use_axes)

    def open(self):
        return self._open(self._filename, 'ab')

    def write_header(self, headers):
        """headers -- string generator"""
        file = self.open()
        for h in headers:
            file.write(h)
            file.write('\n')
        file.close()

    def write(self, lines):
        """lines -- string generator"""
        self._counter += 1
        file = self.open()
        for line in lines:
            file.write(line)
            file.write('\n')
        file.close()

        if self._counter >= self._nfreq:
            self._ifreq += 1
            self._counter = 0
            self.make_filename(self._ifreq)

    def make_filename(self, ifreq, use_axes=False, digit=5):
        # The .gz is not put if AxisWriter is used or compress=no
        ext = '.gz' if (self._use_zip and not(use_axes)) else ''
        num_fmts = '{:0' + str(digit) + '}'
        self._filename = self._base_fn + num_fmts.format(ifreq) + ext


class AxisWriter(Writer):
    """Called to write files if the calculator gives vectors (heat-flux)
    
    Used in FluxWriter.
    For one frame, writes one file per axis.
    
    Keyword arguments:
    filename -- string, will precede the name of written files
    frequency -- integer 
    compresslevel -- integer, ==0 when compress=no in .cfg file
    axes -- list of strings, each the name of an axis
    """

    def __init__(self, filename, frequency, axes, compresslevel=6):
        Writer.__init__(self, filename, frequency,
                        compresslevel, use_axes=True)

        self._axes = axes
        self._dim = len(axes)
        self.make_afnames()

    def open(self):
        # Overwriting open, this way each file is opened
        return map(self._open, self._afnames, ['ab']*self._dim)

    def close(self, file_):
        file_.close()

    def write_header(self, headers):
        files = self.open()
        for h in headers:
            for file_ in files:
                file_.write(h + '\n')
        map(self.close, files)

    def write(self, lines):
        """
        lines: a generator of list of strings.
               Each string is the line for one axis.
        """
        self._counter += 1
        files = self.open()
        for line in lines:
            for i, file_ in enumerate(files):
                file_.write(line[i] + '\n')
                #files[axis].write('\n')
        map(self.close, files)

        if self._counter >= self._nfreq:
            self._ifreq += 1
            self._counter = 0
            self.make_filename(self._ifreq, use_axes=True)
            self.make_afnames()

    def make_afnames(self):
        # Stores the filenames for each axis
        ext = '.gz' if self._use_zip else ''
        self._afnames = [self._filename + axis + ext for axis in self._axes]


"""
MultiXxWriter: choose the options for Current, Flux and NetCDFFluxWriter.
"""

class MultiCurrentWriter:

    def __init__(self, setting, decomp_list, target_anames, group_names):

        self.__atm_writer = CurrentWriter( setting, decomp_list,
                target_anames, revision='atm')
        self.__grp_writer = CurrentWriter( setting, decomp_list,
                group_names, revision='grp')
        self.__inn_writer = None
        self.__out_writer = None

        if setting.curp.decomp_group_current:
            self.__inn_writer = CurrentWriter( setting, decomp_list,
                    group_names, revision='inn')
            self.__out_writer = CurrentWriter( setting, decomp_list,
                    group_names, revision='out')

    def write_header(self):

        self.__atm_writer.write_header()
        self.__grp_writer.write_header()

        if self.__inn_writer:
            self.__inn_writer.write_header()

        if self.__out_writer:
            self.__out_writer.write_header()

    def write(self, istep, key_to_acurs, key_to_gcurs,
            key_to_icurs=None, key_to_ocurs=None):

        self.__atm_writer.write(istep, key_to_acurs)
        self.__grp_writer.write(istep, key_to_gcurs)

        if self.__inn_writer:
            self.__inn_writer.write(istep, key_to_icurs)

        if self.__out_writer:
            self.__out_writer.write(istep, key_to_ocurs)


class MultiFluxWriter:

    def __init__(self, setting, decomp_list,
                 target_anames, group_names, gpair_table=None, axes=None):

        fmt = setting.output.format 
        grain = setting.curp.flux_grain

        if (grain, fmt) == ('atom', 'ascii'):
            self.__atm_writer = FluxWriter( setting, decomp_list,
                    target_anames, axes, revision='atm' )
            self.__grp_writer = None

        elif (grain, fmt) == ('group', 'ascii'):
            self.__atm_writer = None
            self.__grp_writer = FluxWriter( setting, decomp_list,
                    group_names, axes, gpair_table, revision='grp' )

        elif (grain, fmt) == ('both', 'ascii'):
            self.__atm_writer = FluxWriter( setting, decomp_list,
                    target_anames, axes, revision='atm' )
            self.__grp_writer = FluxWriter( setting, decomp_list,
                    group_names, axes, gpair_table, revision='grp' )

        elif (grain, fmt) == ('atom', 'netcdf'):
            self.__atm_writer = NetCDFFluxWriter( setting, decomp_list,
                    target_anames, axes, revision='atm')
            self.__grp_writer = None

        elif (grain, fmt) == ('group', 'netcdf'):
            self.__atm_writer = None
            self.__grp_writer = NetCDFFluxWriter( setting, decomp_list,
                    group_names, gpair_table, axes, revision='grp')

        elif (grain, fmt) == ('both', 'netcdf'):
            self.__atm_writer = NetCDFFluxWriter( setting, decomp_list,
                    target_anames, axes, revision='atm')
            self.__grp_writer = NetCDFFluxWriter( setting, decomp_list,
                    group_names, gpair_table, axes, revision='grp')

        else:
            pass

    # elif method == 'energy-flux' and fmt=='netcdf':
        # obj = NetCDFFluxWriter(setting=setting, decomp_list, target_anames,
                # group_names=group_names, gpair_table=gpair_table)

        
    def write_header(self):
        if self.__atm_writer: self.__atm_writer.write_header()
        if self.__grp_writer: self.__grp_writer.write_header()

    def write(self, istep, flux_atm=None, flux_grp=None):
        if self.__atm_writer: self.__atm_writer.write(istep, flux_atm)
        if self.__grp_writer: self.__grp_writer.write(istep, flux_grp)

"""
Writers below parse the lines for files, given output of a calculator.
Writer or AxisWriter are then used to write the files in case of ascii.
NetCDFFluxWriter does this part by itself.
"""


class CurrentWriter:

    _flag_char = '%'

    def __init__(self, setting, decomp_list, names, title='', revision=''):
        self.__setting = setting
        self.__revision = revision
        self.__names = names
        self.__title = title

        self._time_fmt = '{:>15.3f} [ps]'
        self._data_fmt = ('{{name:>14s}}  '
               '{{xx:{fmt}}} {{xy:{fmt}}} {{xz:{fmt}}} '
               '{{yx:{fmt}}} {{yy:{fmt}}} {{yz:{fmt}}} '
               '{{zx:{fmt}}} {{zy:{fmt}}} {{zz:{fmt}}} ').format(fmt='012.7e')

        # Setup the object using setting
        filename = setting.output.filename[0]
        frequency = setting.output.frequency
        if setting.output.compress:
            compresslevel = 6
        else:
            compresslevel = 0

        prefix, ext = os.path.splitext(filename)
        fmt = '{prefix}_' + self.__revision + '{ext}'
        mod_fn = fmt.format(prefix=prefix, ext=ext)
        self.__key_to_writer = dict(
                total = Writer(mod_fn, frequency, compresslevel) )

        self.__decomp_list = decomp_list
        if setting.output.decomp:
            for key in decomp_list + ['kinetic']:
                key_fn = self.get_key_filename(key, filename)
                self.__key_to_writer[key] = Writer(
                        key_fn, frequency, compresslevel)

    def get_key_filename(self, key, base_filename):
        prefix, ext = os.path.splitext(base_filename)
        fmt = '{prefix}_' + self.__revision + '_{key}{ext}'
        return fmt.format(prefix=prefix, key=key, ext=ext)

    def get_header_format(self, title):
        yield self._flag_char + 'title  ' + title
        yield self._flag_char + 'format ' + 'time ' + self._time_fmt
        yield self._flag_char + 'format ' + 'data ' + self._data_fmt

    def format(self, time, current, names):
        yield self._flag_char + 'time' + self._time_fmt.format(time)
        yield self._flag_char + 'data'
        for name, cur in zip(names, current):
            yield self._data_fmt.format(name=name,
                    xx=cur[0,0], xy=cur[0,1], xz=cur[0,2],
                    yx=cur[1,0], yy=cur[1,1], yz=cur[1,2],
                    zx=cur[2,0], zy=cur[2,1], zz=cur[2,2] )

    def write_header(self):
        # headers
        for key, writer in self.__key_to_writer.items():
            writer.write_header(self.get_header_format(self.__title))

    def write(self, istep, results):

        if self.__setting.output.decomp:
            for key in self.__decomp_list + ['kinetic']:
                lines = self.format(istep, results[key], self.__names)
                self.__key_to_writer[key].write(lines)

        # if self.__setting.output.decomp:
        #     self.
        #     if self.__setting.output.decomp_bonded:
        #         bond_types = ['bond','angle','torsion','improper']
        #     else:
        #         bond_types = ['bonded']
        #         for key in 
        #             lines = self.format(time, results[key], results['names'])
        #             self.__key_to_writer[key].write(lines)

        # else:
        #     for key in ['bonded', 'coulomb', 'vdw', 'vel']:
        #         lines = self.format(time, results[key], results['names'])
        #         self.__key_to_writer[key].write(lines)

        total_current = results['kinetic']
        for key in self.__decomp_list:
            # total_current = total_current + results[key]
            total_current += results[key]

        lines = self.format(istep, total_current, self.__names)
        self.__key_to_writer['total'].write(lines)

    def write_all(self, time_itr, results_itr):
        for time, results in zip(time_itr, results_itr):

            if self.__setting.output.decomp:
                for key in self.__decomp_list + ['kinetic']:
                    lines = self.format(time, results[key], self.__names)
                    self.__key_to_writer[key].write(lines)

            # if self.__setting.output.decomp:
            #     self.
            #     if self.__setting.output.decomp_bonded:
            #         bond_types = ['bond','angle','torsion','improper']
            #     else:
            #         bond_types = ['bonded']
            #         for key in 
            #             lines = self.format(time, results[key], results['names'])
            #             self.__key_to_writer[key].write(lines)

            # else:
            #     for key in ['bonded', 'coulomb', 'vdw', 'vel']:
            #         lines = self.format(time, results[key], results['names'])
            #         self.__key_to_writer[key].write(lines)

            total_current = results['kinetic']
            for key in self.__decomp_list:
                # total_current = total_current + results[key]
                total_current += results[key]

            lines = self.format(time, total_current, names)
            self.__key_to_writer['total'].write(lines)

    def gen_decomp_keys(self, decomp_list):
        # Generate flatten decomposition key list.
        for key in decomp_list:
            if isinstance(key, tuple) or isinstance(key, list):
                yield key[0]
            else:
                yield key


class FluxWriter:

    _flag_char = '%'

    def __init__(self, setting, decomp_list, names, axes,
                 pair_table=None, title='', revision=''):
        self.__setting = setting    # setting object from .cfg file
        self.__revision = revision
        self.__title = title
        self.__names = names
        self.__pair_table = pair_table  # Interaction table
        self.__time_fmt = '{:>15.3f} [ps]'

        if axes is None:
            self.__axes = ['']
        else:
            self.__axes = axes  # List of strings: the name of the axes

        self.__dim = len(self.__axes)

        # Setup the object using setting

        filename = setting.output.filename[0]

        frequency = setting.output.frequency
        if setting.output.compress:
            compresslevel = 6
        else:
            compresslevel = 0
        prefix, ext = os.path.splitext(filename)
        fmt = '{prefix}_' + self.__revision + '{ext}'
        mod_fn = fmt.format(prefix=prefix, ext=ext)
        
        self.__writer = AxisWriter(mod_fn, frequency, self.__axes,
                                   compresslevel)
            
        # Defines format
        # Header columns:
        self.__fmt_fmt = '{:>5} {:>12}'.format('donor', 'acceptor')
        # Values lines format:
        self.__data_fmt = '{:>12s} {:>12s}' 
        # Stores the name of the forces if output has to be decomposed
        self.__decomps = ['total']
        if setting.output.decomp:
            self.__decomps.extend(decomp_list)

        for ptype in self.__decomps:
            self.__fmt_fmt += ' {:>16}'.format(ptype)
        
        self.__num_decomps = len(self.__decomps)
        self.__data_fmt += ' {:>16.8e}' * self.__num_decomps

        # make {gname : igrp} dictionary
        self.__name_to_idx = None
        if pair_table:
            self.__name_to_idx = {
                        name: i for i, name in enumerate(self.__names)
                        }

    def get_header_format(self, title):
        yield self._flag_char + 'title  ' + title
        yield self._flag_char + 'format ' + 'time ' + self.__time_fmt
        yield self._flag_char + 'label ' + self.__fmt_fmt

    def write_header(self):
        self.__writer.write_header(self.get_header_format(self.__title))

    def write(self, istep, key_to_fluxes):
    
        to_fluxes_ar = numpy.array(key_to_fluxes)
        if self.__pair_table is None:
            lines = self.format(istep, key_to_fluxes, self.__names)
        else:
            lines = self.format_by_pairs(
                    istep, key_to_fluxes, self.__names, self.__pair_table)

        self.__writer.write(lines)

    def ok_write(self, i, j):
        if i <= j:
            return True
        else:
            return self.__setting.curp.enable_inverse_pair

    def format(self, time, key_to_fluxes, names):
        """Format each line, returns 2d numpy array"""
        yield [self._flag_char + 'time' + self.__time_fmt.format(time)] * 3
        yield [self._flag_char + 'data'] * 3
        # for name, current in zip(names, current):

        for itar_1, name_i in enumerate(names):
            for jtar_1, name_j in enumerate(names):
                if not self.ok_write(itar_1, jtar_1): continue

                """fluxes = [None for k in self.__decomps]
                for i, pot_type in enumerate(self.__decomps):
                    fluxes[i] = key_to_fluxes[pot_type][itar_1, jtar_1]

                fluxes_ar = numpy.array(fluxes).transpose()
                fluxes_ar = fluxes_ar.reshape(self.__dim, -1)

                line = [self.__data_fmt.format(name_i, name_j, *k)
                        for k in fluxes_ar]
                """

                fluxes = numpy.zeros([self.__num_decomps, self.__dim])
                for i, pot_type in enumerate(self.__decomps):
                    fluxes[i] = key_to_fluxes[pot_type][itar_1, jtar_1]

                line = [self.__data_fmt.format(name_i, name_j, *fluxes[:,i])
                        for i, axis in enumerate(self.__axes)]
                yield line

    def format_by_pairs(self, time, key_to_fluxes, names, pair_table):
        yield [self._flag_char + 'time' + self.__time_fmt.format(time)] * 3
        yield [self._flag_char + 'data'] * 3
        # for name, current in zip(names, current):
        for name_i, names_j in pair_table:
            itar_1 = self.__name_to_idx[name_i]

            for name_j in names_j:
                jtar_1 = self.__name_to_idx[name_j]
                if not self.ok_write(itar_1, jtar_1): continue

                fluxes = numpy.zeros([self.__num_decomps, self.__dim])
                for i, pot_type in enumerate(self.__decomps):
                    fluxes[i] = key_to_fluxes[pot_type][itar_1, jtar_1]
                
                """
                fluxes is either a list of floats or of list.
                Each float corresponding to the contribution of the forces
                (total, vdw, coulomb, ...) to the flux (J).
                Or each list containing floats, coordinates of the flux vector,
                each list element correspond to the contribution of the forces.
                """

                line = [self.__data_fmt.format(name_i, name_j, *fluxes[:,i])
                        for i, axis in enumerate(self.__axes)]
                yield line

    def write_all(self, time_itr, results_itr):
        for time, results in zip(time_itr, results_itr):

            if self.__setting.output.decomp:
                for key in self.__decomp_list + ['kinetic']:
                    lines = self.format(time, results[key], self.__names)
                    self.__key_to_writer[key].write(lines)

            # if self.__setting.output.decomp:
            #     self.
            #     if self.__setting.output.decomp_bonded:
            #         bond_types = ['bond','angle','torsion','improper']
            #     else:
            #         bond_types = ['bonded']
            #         for key in 
            #             lines = self.format(time, results[key], results['names'])
            #             self.__key_to_writer[key].write(lines)

            # else:
            #     for key in ['bonded', 'coulomb', 'vdw', 'vel']:
            #         lines = self.format(time, results[key], results['names'])
            #         self.__key_to_writer[key].write(lines)

            total_current = results['kinetic']
            for key in self.__decomp_list:
                # total_current = total_current + results[key]
                total_current += results[key]

            lines = self.format(time, total_current, self.__names)
            self.__key_to_writer['total'].write(lines)

    def gen_decomp_keys(self, decomp_list):
        """Generate flatten decomposition key list."""
        for key in decomp_list:
            if isinstance(key, tuple) or isinstance(key, list):
                yield key[0]
            else:
                yield key


import netCDF4 as netcdf
class NetCDFFluxWriter:
    """

    For example:

    netcdf file_name {
    dimensions:
            nframe = UNLIMITED ; // (20 currently)
            ngroup_pair = 3 ;
            ncomponent = 9 ;
    variables:
            float time(nframe) ;
                    time:units = "picosecond" ;
            string donors(ngroup_pair) ;
            string acceptors(ngroup_pair) ;
            string components(ncomponent) ;
            float flux(nframe, ngroup_pair, ncomponent) ;
                    flux:units = "kcal/mol/fs" ;

    // global attributes:
                    :title = "" ;
                    :application = "the CURP program" ;
                    :program = "curp" ;
                    :programVersion = "0.7" ;
                    :Convetsions = "CURP" ;
                    :ConvetsionVersion = "0.7" ;
    data:

     time = 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 
        0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19 ;

     donors = "00001_ALA", "00001_ALA", "00002_ALA" ;

     acceptors = "00002_ALA", "00003_ALA", "00003_ALA" ;

     components = "total", "bond", "angle", "torsion", "improper", "coulomb14", 
        "vdw14", "coulomb", "vdw" ;

     flux =
      0.3240475, -0.02477073, 0.01304006, 0.08249603, 0.2595297, -0.06853036, 
        0.01592382, 0.05231608, -0.005957155,
        .
        .
        .,
    }
    """

    version = 0.8

    def __init__(self, setting, decomp_list, names,
                 gpair_table=None, axes=None, title='', revision=''):

        self.__ncfile = None
        self.__setting = setting
        self.__revision = revision
        self.__title = title
        self.__names = names
        self.__pair_table = gpair_table
        self.__fp = setting.output.filename[0]
        self.__axes = axes
        # self.complevel = complevel
        # if setting.output.compress:
            # self.complevel = 1
            # self.use_zlib  = True
        # else:
            # self.complevel = 0
            # self.use_zlib  = False

        self.setup_components(decomp_list)

        self.setup()

    def setup(self):

        fp = self.get_fp()
        title = self.__title

        self.__ngroup_pair = len(self.__names)**2

        ncfile = netcdf.Dataset(fp, clobber=True,
                mode='w', format='NETCDF3_64BIT')

        donors, acceptors = self.setup_gpairs(self.__names, self.__pair_table)

        # set global attributes
        ncfile.title             = title
        ncfile.application       = 'the CURP program'
        ncfile.program           = 'curp'
        ncfile.programVersion    = str(self.version)
        ncfile.Convetsions       = 'CURP'
        ncfile.ConvetsionVersion = str(self.version)

        # create dimensions
        ncfile.createDimension('nframe', None)
        ncfile.createDimension('npair', self.__ngroup_pair)
        ncfile.createDimension('ncomponent' , len(self.__decomps))
        ncfile.createDimension('nchar', 20)
        # file.createDimensions('label', ?)

        # create variables
        # time
        nc_time = ncfile.createVariable('time', 'f4', ('nframe',))
        nc_time.units = 'picosecond'

        # donor and acceptor
        nc_don = ncfile.createVariable('donors',    'c', ('npair', 'nchar'))
        nc_acc = ncfile.createVariable('acceptors', 'c', ('npair', 'nchar'))

        # components
        nc_com = ncfile.createVariable('components','c',
                                      ('ncomponent','nchar'))

        """
        Differences whether the data format is
         - a float (energy-flux) or
         - a vector (heat-flux)
        """
        if self.__axes == None:
            # flux trajectory
            nc_flux = ncfile.createVariable('flux', 'f4',
                    ('nframe', 'npair', 'ncomponent'), )
            self.write = self.energy_write
        else:  
            # create axes dimension
            ncfile.createDimension('naxes', len(self.__axes))
            # flux trajectory
            nc_flux = ncfile.createVariable('flux', 'f4',
                    ('nframe', 'naxes', 'npair', 'ncomponent'), ) 
                    # zlib=self.use_zlib, complevel=self.complevel)
            self.write = self.heat_write

        nc_flux.units = 'kcal/mol/fs'

        # write group pairs
        for ipair_1, (don, acc) in enumerate(zip(donors, acceptors)):
            nc_don[ipair_1] = don.ljust(20)
            nc_acc[ipair_1] = acc.ljust(20)

        # write component names
        for icom_1, com in enumerate(self.__decomps):
            nc_com[icom_1] = com.ljust(20)

        ncfile.sync()
        ncfile.close()
        # self.__ncfile = ncfile

    def setup_gpairs(self, names, pair_table=None):

        name_to_idx = {name:i for i, name in enumerate(names)}
        ntarget = len(names)
        self.__npair_total = ntarget*ntarget

        # make mask indices and name pair list
        mask_indices = []
        donors = []
        acceptors = []
        if pair_table:
            for name_i, names_j in pair_table:
                itar_1 = name_to_idx[name_i]

                for name_j in names_j:
                    jtar_1 = name_to_idx[name_j]
                    if not self.ok_pair(itar_1, jtar_1): continue

                    ipair = ntarget*itar_1 + jtar_1
                    mask_indices += [ipair]
                    donors       += [name_i]
                    acceptors    += [name_j]

        else:
            for name_i in self.__names:
                itar_1 = name_to_idx[name_i]

                for name_j in self.__names:
                    jtar_1 = name_to_idx[name_j]
                    if not self.ok_pair(itar_1, jtar_1): continue

                    ipair = ntarget*itar_1 + jtar_1
                    mask_indices += [ipair]
                    donors       += [name_i]
                    acceptors    += [name_j]

        self.__mask_indices = mask_indices
        self.__ngroup_pair = len(self.__mask_indices)

        return donors, acceptors

    def setup_components(self, decomp_list):

        # Write components
        if self.__setting.output.decomp:
            self.__decomps = ['total'] + decomp_list
        else:
            self.__decomps = ['total']

    def ok_pair(self, i, j):
        if i <= j:
            return True
        else:
            return self.__setting.curp.enable_inverse_pair

    def write_header(self):
        pass

    def get_fp(self):
        prefix, ext = os.path.splitext(self.__fp)
        fmt = '{prefix}_' + self.__revision + '{ext}'
        return fmt.format(prefix=prefix, ext=ext)

    def energy_write(self, istp_1, key_to_fluxes):

        ncfile = self.open()

        # Write time
        dt = 0.01
        ncfile.variables['time'][istp_1] = istp_1*dt

        # Write flux
        nc_flux  = ncfile.variables['flux']

        flux = numpy.array([ key_to_fluxes[pot_type].ravel()
                for pot_type in self.__decomps ]).T

        if self.__mask_indices:
            nc_flux[istp_1] = flux[self.__mask_indices,:]
        else:
            nc_flux[istp_1] = flux
        
        self.close()

    def heat_write(self, istp_1, key_to_fluxes):

        ncfile = self.open()

        # Write time
        dt = 0.01
        ncfile.variables['time'][istp_1] = istp_1*dt

        # Write flux
        nc_flux  = ncfile.variables['flux']

        flux = numpy.array([ key_to_fluxes[pot_type]
                for pot_type in self.__decomps ]).T.reshape(
                                                    (len(self.__axes),
                                                     -1,
                                                     len(self.__decomps))
                                                    )

        if self.__mask_indices:
            nc_flux[istp_1] = flux[:,self.__mask_indices,:]
        else:
            nc_flux[istp_1] = flux
        
        self.close()
    def open(self):
        if self.__ncfile is None:
            ncfile = netcdf.Dataset(self.get_fp(), mode='r+')
            self.__ncfile = ncfile

        return self.__ncfile

    def close(self):
        if self.__ncfile:
            self.__ncfile.sync()
            self.__ncfile.close()
            self.__ncfile = None

