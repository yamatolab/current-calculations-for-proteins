from __future__ import print_function

import os, sys
import math
import numpy
from math import sqrt, pi
from abc import abstractmethod, abstractproperty, ABCMeta

topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)
import utility

curp_dir = os.environ['CURP_HOME']
script_dir = os.path.abspath(os.path.join(curp_dir, 'script'))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)
import pdb_ as PDB

# curp modules
import clog as logger
from exception import CurpException

class NotFoundPropertyError(CurpException): pass
class NumberOfTargetError(CurpException): pass
class VolumeNotDefinedError(CurpException): pass
class TestError(CurpException): pass
class NumberOfGroupError(CurpException): pass

def nproperty(function):
    """Decorator to easy use property.

    Usage:

    class A(object):
        def __init__(self):
            self.__x = None

        @nproperty
        def x():

            test_val_list = [] # use in get and set

            def get(self):
                print('Get property x.')
                return self.__x

            def set(self, x):
                print('Set property x: ', x, '.')
                self.__x = x

            def del_(self, x):
                del self.__x

    a = A()
    a.x = 5
    print a.x
    => 5
    getattr(a, 'x') => 5
    setattr(a, 'x', 9)
    print a.x
    => 9
    """

    func_locals = {'doc':function.__doc__}
    def probe_func(frame, event, arg):
        if event == 'return':
            locals = frame.f_locals

            prop_dict = {}
            for name, method in locals.items():
                if name.startswith('get'):
                    prop_dict['fget'] = method
                    get_method = method
                elif name.startswith('set'):
                    prop_dict['fset'] = method
                elif name.startswith('del'):
                    prop_dict['fdel'] = method
                # elif name.startswith('doc'):
                #     prop_dict['doc']  = method
                else:
                    pass

            func_locals.update( prop_dict )
            sys.settrace(None)
        return probe_func
    sys.settrace(probe_func)
    function()
    return property(**func_locals)

class UndefinedPropertyError(CurpException): pass
class VolumeSetting:

    def __init__(self, topology=None, setting=None,
            traj_parser=None, target_atoms=None, gname_iatoms_pairs=None):
        self.__tpl = topology
        self.__setting = setting
        self.__traj_parser = traj_parser
        self.__target_atoms = target_atoms
        self.__gname_iatoms_pairs = gname_iatoms_pairs

        self.__props = dict(
                natom               = 0  , 
                names               = [] , 
                elems               = [] , 
                resnames            = [] , 
                gname_iatoms_pairs  = [] , 
                target_atoms        = [] , 
                output_volume_file  = '' , 
                output_gvolume_file = '' , 

                vdw_radii = [] , 

                smve_rmax      = 2.5  , 
                smve_dr        = 0.01 , 
                smve_interval  = 1    , 
                smve_increment = 5    , 

                atomic_trajectory_file = '' , 
                group_trajectory_file  = '' , 

                voronoi_cutoff       = 6.0   , 
                voronoi_no_hydrogen  = False , 
                voronoi_solvation    = 'none' , 
                voronoi_probe_length = 2.4   , 
                voronoi_output_solvation_file = '',

                group_method = 'none'
        )


    def __str__(self):
        return str(self.__props)

    @nproperty
    def traj_parser():

        def get(self):
            return self.__traj_parser

        def set(self, value):
            self.__traj_parser = value

    @nproperty
    def natom():

        def get(self):
            if self.__tpl:
                return self.__tpl.get_natom()
            else:
                return self.__props['natom']

        def set(self, value):
            if self.__tpl is None:
                self.__props['natom'] = value

    @nproperty
    def names():

        def get(self):
            if self.__tpl:
                return self.__tpl.get_atom_info()['names']
            else:
                return self.__names

        def set(self, value):
            if self.__tpl is None:
                self.__props['names'] = value

    @nproperty
    def elems():

        def get(self):
            if self.__tpl:
                return self.__tpl.get_atom_info()['elems']
            else:
                return self.__elems

        def set(self, value):
            if self.__tpl is None:
                self.__props['elems'] = value

    @nproperty
    def resnames():

        def get(self):
            if self.__tpl:
                return self.__tpl.get_residue_info()['names']
            else:
                return self.__resnames

        def set(self, value):
            if self.__tpl is None:
                self.__props['resnames'] = value

    @nproperty
    def rids():

        def get(self):
            if self.__tpl:
                return self.__tpl.get_residue_info()['ids']
            else:
                return self.__rids

        def set(self, value):
            if self.__tpl is None:
                self.__props['rids'] = value

    @nproperty
    def target_atoms():

        def get(self):
            if self.__target_atoms:
                return self.__target_atoms
            else:
                return self.__props['target_atoms']

        def set(self, value):
            if self.__target_atoms is None:
                self.__props['target_atoms'] = value

    @nproperty
    def gname_iatoms_pairs():

        def get(self):
            if self.__gname_iatoms_pairs:
                return self.__gname_iatoms_pairs
            else:
                return self.__props['gname_iatoms_pairs']

        def set(self, value):
            if self.__gname_iatoms_pairs is None:
                self.__props['gname_iatoms_pairs'] = value

    @nproperty
    def output_volume_file():

        def get(self):
            if self.__setting:
                return self.__setting.volume.output_volume_file[0]
            else:
                return self.__props['output_volume_file']

        def set(self, value):
            if self.__setting is None:
                self.__props['output_volume_file'] = value

    @nproperty
    def output_gvolume_file():

        def get(self):
            if self.__setting:
                return self.__setting.volume.output_gvolume_file[0]
            else:
                return self.__props['output_gvolume_file']

        def set(self, value):
            if self.__setting is None:
                self.__props['output_gvolume_file'] = value

    @nproperty
    def vdw_radii():

        def get(self):
            if self.__tpl:
                return self.__tpl.get_atom_info()['vdw_radii']
            else:
                return self.__vdw_radii

        def set(self, value):
            if self.__tpl is None:
                self.__props['vdw_radii'] = value

    @nproperty
    def smve_rmax():

        def get(self):
            if self.__setting:
                return self.__setting.volume.smve_rmax
            else:
                return self.__props['smve_rmax']

        def set(self, value):
            if self.__setting is None:
                self.__props['smve_rmax'] = value

    @nproperty
    def smve_dr():

        def get(self):
            if self.__setting:
                return self.__setting.volume.smve_dr
            else:
                return self.__props['smve_dr']

        def set(self, value):
            if self.__setting is None:
                self.__props['smve_dr'] = value

    @nproperty
    def smve_interval():

        def get(self):
            if self.__setting:
                return self.__setting.volume.smve_interval
            else:
                return self.__props['smve_interval']

        def set(self, value):
            if self.__setting is None:
                self.__props['smve_interval'] = value

    @nproperty
    def smve_increment():

        def get(self):
            if self.__setting:
                return self.__setting.volume.smve_increment
            else:
                return self.__props['smve_increment']

        def set(self, value):
            if self.__setting is None:
                self.__props['smve_increment'] = value

    @nproperty
    def atomic_trajectory_file():

        def get(self):
            if self.__setting:
                return self.__setting.volume.atomic_trajectory_file[0]
            else:
                return self.__props['atomic_trajectory_file']

        def set(self, value):
            if self.__setting is None:
                self.__props['atomic_trajectory_file'] = value

    @nproperty
    def group_trajectory_file():

        def get(self):
            if self.__setting:
                return self.__setting.volume.group_trajectory_file[0]
            else:
                return self.__props['group_trajectory_file']

        def set(self, value):
            if self.__setting is None:
                self.__props['group_trajectory_file'] = value

    @nproperty
    def voronoi_no_hydrogen():

        def get(self):
            if self.__setting:
                return self.__setting.volume.voronoi_no_hydrogen
            else:
                return self.__props['voronoi_no_hydrogen']

        def set(self, value):
            if self.__setting is None:
                self.__props['voronoi_no_hydrogen'] = value

    @nproperty
    def voronoi_cutoff():

        def get(self):
            if self.__setting:
                return self.__setting.volume.voronoi_cutoff
            else:
                return self.__props['voronoi_cutoff']

        def set(self, value):
            if self.__setting is None:
                self.__props['voronoi_cutoff'] = value

    @nproperty
    def voronoi_solvation():

        def get(self):
            if self.__setting:
                return self.__setting.volume.voronoi_solvation
            else:
                return self.__props['voronoi_solvation']

        def set(self, value):
            if self.__setting is None:
                self.__props['voronoi_solvation'] = value

    @nproperty
    def voronoi_probe_length():

        def get(self):
            if self.__setting:
                return self.__setting.volume.voronoi_probe_length
            else:
                return self.__props['voronoi_probe_length']

        def set(self, value):
            if self.__setting is None:
                self.__props['voronoi_probe_length'] = value

    @nproperty
    def voronoi_output_solvation_file():

        def get(self):
            if self.__setting:
                return self.__setting.volume.voronoi_output_solvation_file[0]
            else:
                return self.__props['voronoi_output_solvation_file']

        def set(self, value):
            if self.__setting is None:
                self.__props['voronoi_output_solvation_file'] = value

    @nproperty
    def group_method():

        def get(self):
            if self.__setting:
                return self.__setting.curp.group_method
            else:
                return self.__props['group_method']

        def set(self, value):
            if self.__setting is None:
                self.__props['group_method'] = value

class VolumeCalculatorBase:

    __metaclass__ = ABCMeta

    def __init__(self, volume_setting):
        self._setting = volume_setting
        self.__iatm_to_itars = None
        self.__do_write_volume = False
        self.__do_write_gvolume = False

        # prepare to write out volumes and group volumes data
        self.__count = 0
        if self._setting.output_volume_file != '':
            self.__do_write_volume = True

        if self._setting.output_gvolume_file != '':
            self.__do_write_gvolume = True

        self.__gnames = [ gname for gname, iatoms
                in self._setting.gname_iatoms_pairs ]

        self.prepare()

    @abstractmethod
    def prepare(self):
        return

    @abstractmethod
    def cal_volume(self, crd):
        return

    def get_volume(self, crd=[]):
        self.__count += 1
        self.__volumes = self.cal_volume(crd)

        # write volumes
        self.write_atomic_volume(self.__volumes)

        return self.__volumes

    def get_gvolume(self, crd=[]):
        """Sum up the volume in terms with given group and get its volume."""

        if self.__volumes is None:
            return None

        else:
            gvolumes = []
            iatm_to_itars = self.get_iatm_to_itars()

            for igrp_1, (gname, atoms) in enumerate(
                    self._setting.gname_iatoms_pairs):

                gvol = 0.0
                
                for iatm in atoms:
                    itar = iatm_to_itars[iatm-1]
                    if itar == 0: continue

                    gvol += self.__volumes[itar-1]

                gvolumes.append( gvol )

            # write volumes
            self.write_group_volume(gvolumes)

            return gvolumes

    def get_iatm_to_itars(self):
        if self.__iatm_to_itars is None:
            iatm_to_itars = numpy.zeros( [self._setting.natom], numpy.int)

            for itar_1, iatm in enumerate(self._setting.target_atoms):
                iatm_to_itars[iatm-1] = itar_1 + 1

            self.__iatm_to_itars = iatm_to_itars

        return self.__iatm_to_itars

    def get_itar_by_iatm(self, iatm=0):
        iatm_to_itars = self.get_iatm_to_itars()
        natom = len(iatm_to_itars)
        if iatm <= 0:
            return 0
        elif iatm > natom:
            return 0
        else:
            return iatm_to_itars[iatm-1]

    def write_atomic_volume(self, volumes):
        if not self.__do_write_volume: return

        filename = self._setting.output_volume_file
        names = self._setting.names

        fmt = '{id:>5d} {name:<10s}   {volume:>12.5f}'
        lines = []

        # for system
        file = open(filename, 'a')
        file.write('#BEGIN step ' + str(self.__count) + '\n')

        for iatm, volume in zip(self._setting.target_atoms, volumes):
            name = names[iatm-1]
            line = fmt.format(id=iatm, name=name, volume=volume)
            file.write(line + '\n')

        file.write('#END\n')
        file.close()

    def write_group_volume(self, gvolumes):
        if not self.__do_write_gvolume: return

        filename = self._setting.output_gvolume_file
        gnames = self.__gnames

        fmt = '{id:>5d} {name:<10s}   {volume:>12.5f}'
        lines = []

        # for system
        file = open(filename, 'a')
        file.write('#BEGIN step ' + str(self.__count) + '\n')
        logger.debug(gnames)
        logger.debug(gvolumes)

        for i_1, (name, volume) in enumerate( zip(gnames, gvolumes) ):
            line = fmt.format(id=i_1+1, name=name, volume=volume)
            file.write(line + '\n')

        file.write('#END\n')
        file.close()


class VDWVolumeCalculator(VolumeCalculatorBase):

    def __init__(self, volume_setting):
        self.__volumes = None
        VolumeCalculatorBase.__init__(self, volume_setting)

    def prepare(self):
        self.__target_atoms = self._setting.target_atoms
        self.__vdw_radii = self._setting.vdw_radii
        self.__volumes = self.cal_volume()

    def cal_volume(self, crd=[]):
        """Calculate the atom volumes from atom radii."""

        if self.__volumes is None:
            volumes = []
            for iatm in self.__target_atoms:
                radius = self.__vdw_radii[iatm-1]
                volumes.append( 4*pi/3 * radius**3 )
            self.__volumes = numpy.array(volumes)

        return self.__volumes

class SMVEVolumeCalculator(VolumeCalculatorBase):

    """
    Method to calculate the volume introduced by Srolovitz, Maeda, 
    Vitek and Egami.
    """

    def __init__(self, volume_setting):
        VolumeCalculatorBase.__init__(self, volume_setting)

    def prepare(self):
        self.__well2s = numpy.array(list(self.gen_wells()))**2

    def gen_wells(self):
        """Generate the well from radial distribution functions."""

        # get parameter from setting
        traj_parser = self._setting.traj_parser
        rmax     = self._setting.smve_rmax
        dr       = self._setting.smve_dr
        interval = self._setting.smve_interval
        incr     = self._setting.smve_increment

        # calculate radial distribution functions
        from calrdf import average_rdf
        rdfs = average_rdf(traj_parser, rmax=rmax, dr=dr, interval=interval,
                average=True, per_area=True) 

        # search 1st wells
        from search_well import gen_searched_wells
        for well_info in gen_searched_wells(rdfs, dr, increment=incr):
            tmp1, tmp2 ,well, tmp4 = well_info
            yield well

    def get_volume(self, crd):
        radii, volumes = self.get_volume_fort(crd)

        # for iatm_1, (radius, volume) in enumerate(zip(radii, volumes)):
        #     logger.info(iatm_1+1, radius, volume)

        return volumes

    def get_volume_py(self, crd):
        """Calculate and get the volume by python code."""
        radii = numpy.array( list(self.gen_radius(crd)) )
        volumes = 4.0*math.pi*radii**3 / 3.0
        return radii, volumes

    def gen_radius(self, crd):
        """Calculate a radius and volume on each atoms by wells."""
        # wells2 is powered by its self(wells).

        for iatm_1 in range(self._setting.natom):
            well2 = self.__well2s[iatm_1]
            num = 0.0 # numerator
            den = 0.0 # denominator

            for jatm_1 in range(self._setting.natom.natom):
                if iatm_1 == jatm_1: continue
                r_ij = crd[iatm_1] - crd[jatm_1]
                l_ij2 = numpy.dot(r_ij, r_ij)

                if l_ij2 > well2: continue

                l_ij = math.sqrt(l_ij2)

                num += 1.0/l_ij
                den += 1.0/l_ij2

            yield 0.5 * num / den

    def get_volume_fort(self, crd):
        """Calculate and get the volume by fortran code."""
        from lib_calvolume import calvolume
        radii, volumes = calvolume(crd, self.__well2s)
        return radii, volumes

import lib_voronoi as lib_voro
import lib_sink
class VoronoiVolumeCalculatorBase(VolumeCalculatorBase):
    
    """
    This class is the base class to use o calculate the volume obtained
    from voronoi polyhedra.
    Please beware because this method is time-consuming routine.
    """

    def __init__(self, vsetting):
        VolumeCalculatorBase.__init__(self, vsetting)

    @abstractmethod
    def preprocess_crd(self, crd):
        """Preprocess the coordinate by any way."""
        return

    @abstractmethod
    def get_enable_atoms(self, crd):
        """Get the boolean array that atoms are enable for voronoi method."""
        return

    def cal_volume(self, crd):
        """Calculate and get the volume by fortran code."""
        self.__radii, self.__volumes = self.cal_voronoi(crd)
        # radii, volumes = self.cal_voronoi_with_debug(crd)
        return self.__volumes

    def init_voronoi(self, crd):
        """Initialize voronoi calculation for one snapshop."""

        # preprocess
        new_crd = self.preprocess_crd(crd)

        # get enable atoms for current coordinate
        self.__is_enable_atoms = self.get_enable_atoms(new_crd)

        # get cutoff
        cutoff = self._setting.voronoi_cutoff

        # initialize the voronoi_one module
        lib_voro.voronoi_one.initialize(self.__is_enable_atoms, cutoff, 30.0)

        return new_crd, self.__is_enable_atoms

    def get_gvolume(self, crd=[]):
        nohyd = self._setting.voronoi_no_hydrogen
        if self._setting.group_method == 'united' and nohyd:
            gvolumes = []

            # for itar_1, iatm in enumerate(self._setting.target_atoms):
            #     if not self.__is_enable_atoms[iatm-1]: continue

            #     gvolumes.append( self.__volumes[itar_1] )

            for gname, iatoms in self._setting.gname_iatoms_pairs:
                gvol = None
                for iatm in iatoms:
                    if self.__is_enable_atoms[iatm-1]:
                        itar = self.get_itar_by_iatm(iatm)
                        gvol = self.__volumes[itar-1]
                        break

                if gvol is None:
                    msg = ('At gname = {} volume cannot be defined. '
                           'Please set the difinition of target_atoms or group '
                           'to the appropriate value.')
                    raise VolumeNotDefinedError( msg.format(gname) )
                else:
                    gvolumes.append( gvol )

            # write volumes
            self.write_group_volume(gvolumes)

            return gvolumes

        else:
            return VolumeCalculatorBase.get_gvolume(self, crd)
    
    def cal_voronoi(self, crd):
        """Calculate and get the volume by fortran code."""
        # print(1)
        try:
            # print(2)
            crd, enable_atoms = self.init_voronoi(crd)
            # print(3)
        except ValueError:
            # print(4)
            msg = ("When a solvated system is handled, "
                   "the keyword voronoi_solvation should be 'none'.")
            # print(5, msg)
            raise CurpException(msg)
            # print(6)
        # print(7)

        radii, volumes = lib_voro.cal_voronoi_all(crd,
                self._setting.target_atoms, enable_atoms)
        return radii, volumes

    def cal_voronoi_with_debug(self, crd):

        crd, enable_atoms = self.init_voronoi(crd)

        # For the atomic calculation, the case of no hydrogen gives
        #the default value of the hydrogen volume.
        hyd_volume = 8.0

        # give the coordinate to the voronoi_one module
        lib_voro.voronoi_one.crd = crd

        volumes   = []
        nnabs     = []
        nabs_list = []

        ntar = len(self._setting.target_atoms)

        for iatm in self._setting.target_atoms:

            if self.__is_enable_atoms[iatm-1]:
                volume = lib_voro.cal_voronoi_one(iatm)
                self.write_informations(lib_voro.voronoi_one)

                nnabs.append(lib_voro.voronoi_one.nnab)
                nabs_list.append(lib_voro.voronoi_one.nab_list)

            else:
                volume = hyd_volume
                nnabs.append(None)
                nabs_list.append(None)

            volumes.append( volume )

        self.write_summary(nnabs, nabs_list)

        radii = lib_voro.cal_radii_sphere(numpy.array(volumes))
        return radii, volumes

    def write_informations(self, mod):
        """Write out the information of one atom 
        obtained by Voronoi analysis."""

        # write out results
        logger.debug(' NUMBER OF NEIGHBOURS', mod.nface)
        logger.debug(' NEIGHBOUR LIST')

        logger.debug(' {:5s}   {:5s}'.format('ATOM', 'FACE'))
        logger.debug(' {:5s}   {:5s}   {:>36s}   {:>12s}'.format(
                'INDEX', 'EDGES', 'RELATIVE POSITION', 'DISTANCE'))

        can_fmt = " {:5s}   {:5s}   {:12.5f}{:12.5f}{:12.5f}   {:12.5f}"
        for ican_1 in range(mod.ncan):
            if mod.edges[ican_1] == 0: continue

            logger.debug(can_fmt.format(
                mod.can_to_iatm[ican_1], mod.edges[ican_1],
                mod.can_crd[ican_1, 0], mod.can_crd[ican_1, 1],
                mod.can_crd[ican_1, 2], sqrt(mod.can_len2[ican_1])
            ))

        logger.debug(' NUMBER OF EDGES ', mod.nedge)
        logger.debug(' NUMBER OF VERTICS ', mod.nver)
        logger.debug(' VERTEX LIST ')
        logger.debug(' {:>18s}   {:>36s}'.
                format('INDICES', 'RELATIVE POSITION'))

        for iver_1 in range(mod.nver):

            ican = mod.ver_to_3cans[iver_1, 0]
            jcan = mod.ver_to_3cans[iver_1, 1]
            kcan = mod.ver_to_3cans[iver_1, 2]

            # logger.debug(' {:6d}{:6d}{:6d}   {:12.5f}{12.5f}{12.5f})'.format(
            #     mod.can_to_iatm[ican-1], mod.can_to_iatm[jcan-1],
            #     mod.can_to_iatm[kcan-1], mod.ver_crd[iver_1, 0],
            #     mod.ver_crd[iver_1, 1],  mod.ver_crd[iver_1, 2]
            # ))
            logger.debug(' {:6}{:6}{:6}   {:12.5f}{:12.5f}{:12.5f}'.format(
                mod.can_to_iatm[ican-1], mod.can_to_iatm[jcan-1],
                mod.can_to_iatm[kcan-1], mod.ver_crd[iver_1, 0],
                mod.ver_crd[iver_1, 1],  mod.ver_crd[iver_1, 2]
            ))

    def write_summary(self, nnabs, nabs_list):
        """Write out the summary of Voronoi calculation."""

        # ** summary **
        logger.debug(' FINAL SUMMARY')
        logger.debug(' {:5s}   {:5s}   {}'.format(
            'INDEX', 'NABS', '... NEIGHBOURS INDICES' ))

        ncoord = 0
        for itar_1, (iatm, nnab, nab_list) in enumerate(
                zip(self._setting.target_atoms, nnabs, nabs_list)):

            if nnab is None:
                logger.debug( ' {:5} is not calculated.'.format(iatm) )
                continue

            ncoord = ncoord + nnab

            logger.debug( ' {:5}   {:5}   {}'.
                    format( iatm, nnab, nab_list[:nnab]))

            # ** check that if jatm is a neighbour of iatm **
            # ** then iatm is also a neighbour of jatm     **

            # for inab_1 in zip(self.__target_atoms, ange(nnabs[iatm-1]):

            # for jatm, nnab2 in zip(self.__target_atoms, ange(nnabs[iatm-1]):

            for jatm in nab_list[:nnab]:
                jtar = self.get_itar_by_iatm(jatm)
                if jtar == 0: continue
                found = False

                nab_j_list = nabs_list[jtar-1]
                nnab_j     = nnabs[jtar-1]

                for katm in nab_j_list[:nnab_j]:
                    # exit if there is iatm in the jatm-list at least.
                    if found: break

                    found = iatm == katm

                else:
                    logger.debug(' {:5} is not a neighbour of {:5}'.format(
                        iatm, jatm) )

            # for jatm, nnab2 in zip(self.__target_atoms, ange(nnabs[iatm-1]):
            #     jatm = nab_list[inab_1, iatm_1]
            #     found = False

            #     for jnab_1 in range(nnab[jatm-1]):
            #         # exit if there is iatm in the jatm-list at least.
            #         if found: break

            #         found = iatm == nab_list[jnab_1, jatm-1]

            #     else:
            #         logger.debug(' {:5d} is not a neighbour of {:5d}'.format(
            #             iatm, jatm) )


        coord = float(ncoord/float(self._setting.natom))
        logger.debug(' average coordination number = {:10.5f}'.format(coord))

class VoronoiVolumeCalculator(VoronoiVolumeCalculatorBase):
    
    """
    This class will be used if voronoi_solvation is 'none'.
    """

    def __init__(self, vsetting):
        self.__is_enable_atoms = None
        VoronoiVolumeCalculatorBase.__init__(self, vsetting)

    def prepare(self):
        self.__is_enable_atoms = self.get_enable_atoms()

    def preprocess_crd(self, crd):
        return crd

    def get_enable_atoms(self, crd=[]):
        if self.__is_enable_atoms is None:
            if self._setting.voronoi_no_hydrogen:
                is_enable_atoms = numpy.array(
                        [ elem!='H' for elem in self._setting.elems ] )

            else:
                is_enable_atoms = numpy.ones(
                        [self._setting.natom], numpy.bool )

            self.__is_enable_atoms = is_enable_atoms

        return self.__is_enable_atoms

class VoronoiVolumeCalculatorWithSolvation(VoronoiVolumeCalculatorBase):
    
    """
    This class will be used if voronoi_solvation is not 'none'.
    """

    def __init__(self, vsetting):
        self.__water_crd = None
        self.__is_system_atoms = None
        self.__is_enable_system_atoms = None
        self.__itraj = 0
        self.__count = 0
        self.__do_write_crd = False
        self.__fn_base = None
        self.__fn_ext = None
        VoronoiVolumeCalculatorBase.__init__(self, vsetting)

    def prepare(self):

        # parameter
        nohyd  = self._setting.voronoi_no_hydrogen
        stype  = self._setting.voronoi_solvation

        # load the water molucules data from pdb file
        mod_dir = os.path.dirname(__file__)
        type_to_fns = {
                'RANDOM20' : os.path.join(mod_dir, './random20.pdb.gz'),
                'NPT20'    : os.path.join(mod_dir, './npt20.pdb.gz'   ),
                'NPT25'    : os.path.join(mod_dir, './npt25.pdb.gz'   ),
                'NPT30'    : os.path.join(mod_dir, './npt30.pdb.gz'   ),
        }
        self.__water_crd = self.parse_water_pdb( type_to_fns[stype] )

        # make is_enable_system_atoms
        if nohyd:
            self.__is_enable_system_atoms = numpy.array(
                    [ elem!='H' for elem in self._setting.elems ] )

        else:
            self.__is_enable_system_atoms = numpy.ones(
                    [self._setting.natom], numpy.bool)

        # output solvation coordinate
        filename = self._setting.voronoi_output_solvation_file
        if filename != '':
            self.__do_write_crd = True
            self.__fn_base, self.__fn_ext = os.path.splitext(filename)

    def preprocess_crd(self, crd):
        new_crd = self.solvate_crd(crd)
        if self.__do_write_crd:
            self.write_crd(new_crd)
        return new_crd

    def get_enable_atoms(self, crd):
        if self._setting.voronoi_no_hydrogen:
            is_enable_atoms = numpy.ones( [len(crd)], numpy.bool)
            natom = len(self.__is_enable_system_atoms)
            is_enable_atoms[:natom] = self.__is_enable_system_atoms

        else:
            is_enable_atoms = numpy.ones( [len(crd)], numpy.bool)

        return is_enable_atoms

    def parse_water_pdb(self, filename):
        """Parse the water solvent pdb and get the coordinates and
        elements from it."""

        atoms = list(PDB.System(filename))

        water_crd   = []
        if self._setting.voronoi_no_hydrogen:
            for atom in atoms:
                if atom.element == 'H': continue
                water_crd.append( atom.pos )

        else:
            for atom in atoms:
                water_crd.append( atom.pos )

        return numpy.array( water_crd )

    def solvate_crd(self, crd):
        """Return the coordinate that sunk the system into watar molecules."""
        probe = self._setting.voronoi_probe_length
        if self._setting.voronoi_no_hydrogen:
            solv_crd, natom_solv = lib_sink.sink_nohyd(
                    crd, self.__water_crd, probe )
        else:
            solv_crd, natom_solv = lib_sink.sink( crd, self.__water_crd, probe )
        return solv_crd[:natom_solv]

    def write_crd(self, crd):

        line_template = (
                'ATOM  {id:>5} {name:<4} {rname:>3} {cname:>1}{rid:>4}'
            '    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfact:6.2f}'
            '{other:>4}       {elem:<2} '
        )

        names = self._setting.names
        resnames = self._setting.resnames
        rids = self._setting.rids
        elems = self._setting.elems

        occ   = 1.0
        bfact = 0.0
        cname = ''
        other = ''

        pdblines = []

        # for system
        for iatm_1, (name, rname, rid, pos, elem) in enumerate(
                        zip(names, resnames, rids, crd, elems) ):

            pdbline = line_template.format(
                id=iatm_1+1, name=name, rname=rname, rid=rid, cname=cname,
                x=pos[0], y=pos[1], z=pos[2], occ=occ, bfact=bfact,
                other=other, elem=elem )

            pdblines.append( pdbline )

        pdblines.append( 'TER')
        # for solvent
        iatm_end = iatm_1+1
        rid_end = rid

        if self._setting.voronoi_no_hydrogen:
            natom_sol = len(crd) - len(names) 
            nres_sol = natom_sol

            for index in range(nres_sol):
                iatm = iatm_end + index + 1
                rid  = rid_end + index + 1

                pos = crd[iatm-1]
                pdbline = line_template.format(
                    id=iatm, name='O', rname='WAT', rid=rid, cname=cname,
                    x=pos[0], y=pos[1], z=pos[2], occ=occ, bfact=bfact,
                    other=other, elem='O' )
                pdblines.append( pdbline )
                pdblines.append( 'TER')

        else:
            natom_sol = len(crd) - len(names) 
            nres_sol = natom_sol/3

            for index in range(nres_sol):
                iatm = iatm_end + 3*index + 1
                rid  = rid_end + index + 1

                pos = crd[iatm-1]
                pdbline = line_template.format(
                    id=iatm, name='O', rname='WAT', rid=rid, cname=cname,
                    x=pos[0], y=pos[1], z=pos[2], occ=occ, bfact=bfact,
                    other=other, elem='O' )
                pdblines.append( pdbline )

                pos = crd[iatm-1+1]
                pdbline = line_template.format(
                    id=iatm+1, name='H1', rname='WAT', rid=rid, cname=cname,
                    x=pos[0], y=pos[1], z=pos[2], occ=occ, bfact=bfact,
                    other=other, elem='H' )
                pdblines.append( pdbline )

                pos = crd[iatm-1+2]
                pdbline = line_template.format(
                    id=iatm+2, name='H2', rname='WAT', rid=rid, cname=cname,
                    x=pos[0], y=pos[1], z=pos[2], occ=occ, bfact=bfact,
                    other=other, elem='H' )
                pdblines.append( pdbline )

                pdblines.append( 'TER')
    
        self.__count += 1
        filename = '{}{:05}{}'.format(
                self.__fn_base, self.__count, self.__fn_ext )

        file = open(filename, 'w')
        for line in pdblines:
            file.write(line + '\n')
        file.close()


class Volume1(VolumeCalculatorBase):
    
    """
    Method to use 1 as volume values for the group calculations.
    """

    def __init__(self, volume_setting):
        VolumeCalculatorBase.__init__(self, volume_setting)

    def prepare(self):
        self.__volumes  = numpy.ones( [len(self._setting.target_atoms)] )
        self.__gvolumes = numpy.ones( [len(self._setting.gname_iatoms_pairs)] )

    def cal_volume(self, crd=[]):
        """Return 1.0."""
        return self.__volumes

    def get_gvolume(self, crd=[]):
        """Return 1.0."""
        self.write_group_volume(self.__gvolumes)
        return self.__gvolumes

class OuterVolumeFetcher(VolumeCalculatorBase):
    
    """
    A class that performed the method to get outer volume trajectory file
    already calculeted by any methods.
    """

    end_flag = '#END'
    comment  = '#'

    def __init__(self, volume_setting):
        traj_fn = volume_setting.atomic_trajectory_file
        self.__file = open(traj_fn, 'r')

        grp_traj_fn = volume_setting.group_trajectory_file
        if grp_traj_fn == '':
            self.__gfile = None
        else:
            self.__gfile = open(grp_traj_fn, 'r')

        VolumeCalculatorBase.__init__(self, volume_setting)

    def prepare(self):
        self.__ngrp = len(self._setting.gname_iatoms_pairs)
        self.__ntar = len(self._setting.target_atoms)

    def cal_volume(self, crd):
        try:
            ids, names, vols = self.parse_atomic().next()

            if len(vols) != self.__ntar:
                msg = 'The number of target atoms: {}, but in volume trajectory: {}'
                raise NumberOfTargetError(
                        msg.format(self.__ntar, len(vols)) )

            return vols
        except StopIteration:
            return None

    def get_gvolume(self, crd=[]):
        if self.__gfile is None:
            return VolumeCalculatorBase.get_gvolume(self, crd)
        else:
            try:
                ids, names, vols = self.parse_group().next()

                if len(vols) != self.__ngrp:
                    msg = 'The number of groups: {}, but in volume trajectory: {}'
                    raise NumberOfTargetError(
                            msg.format(self.__ngrp, len(vols)) )

                self.write_group_volume(vols)

                return vols
            except StopIteration:
                return None

    def parse_atomic(self):
        return self.parse_trajectory(self.__file)

    def parse_group(self):
        return self.parse_trajectory(self.__gfile)

    def parse_trajectory(self, file):

        lines = []
        for line in file:
            if line.startswith(self.end_flag):
                yield self.parse_snapshop(lines)
                lines = []
            elif line.startswith('#'):
                pass
            else:
                lines.append(line)

        self.__file.close()

    def parse_snapshop(self, lines, other=False):
        natom = len(lines)
        ids = numpy.zeros((natom),dtype=numpy.int)
        names = numpy.zeros((natom), dtype=numpy.str)
        vols = numpy.zeros((natom))

        for iatm, line in enumerate(lines):
            id, name, vol = line.split()
            ids[iatm] = id
            names[iatm] = name
            vols[iatm] = vol

        return ids, names, vols

        for i, atom in enumerate(self.__system.atoms):

            cname = atom.chain.name if atom.chain else ' '
            if atom.chain:
                if not first:
                    if cname != prev_cname:
                        pdblines.append('TER')
                        prev_cname = cname
                else:
                    first = False
            # else:
            #     if rid != prev_rid:
            #         pdblines.append('TER')
            #         prev_rid = rid

            # id = atom.id
            id = i+1
            rid   = atom.residue.id
            rname = atom.residue.name
            x, y, z = atom.pos
            prev_cname = cname
            # prev_rid = rid

            # format the atom name
            aname = atom.name if atom.name[0].isdigit() else ' ' + atom.name

            line_template = (
                'ATOM  {id:>5} {name:<4} {rname:>3} {cname:>1}{rid:>4}'
                '    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfact:6.2f}'
                '{other:>4}       {elem:<2} '
            )
            pdbline = line_template.format(
                id=id, name=aname, rname=rname, rid=rid, cname=cname,
                x=x, y=y, z=z, occ=atom.occ, bfact=atom.bfact,
                other=atom.other, elem=atom.elem )

            if term_resname:
                if term_resname != rname:
                    pdblines.append('TER')
                    term_resname = ""

            pdblines.append(pdbline)

            if atom.name == 'OXT':
                term_resname = rname
                # pdblines.append('TER')

        else:
            if pdblines != []:
                pdblines.append('TER')
            pass

        return '\n'.join(pdblines)


def pdb_test():

    def save_merged_pdb(pdb_fn, bfactors):

        atoms = list(PDB.System(pdb_fn))
        for atom, bfact in zip(atoms, bfactors):
            atom.bfactor = bfact

        return atoms

    def save_merged_pdb_three(pdb_fn, bfactors):

        atoms = System(PDB.System(pdb_fn))
        for atom, bfact in zip(atoms, bfactors):
            if bfact < 1.0:
                bfact = 0.0
            elif 1.0 <= bfact <= 10.0:
                bfact = 10.0
            elif 10.0 < bfact:
                bfact = 50.0
            else:
                pass

            atom.bfactor = bfact

        return atoms

    def save_pdb(atoms, pdb_fn):

        with  open(pdb_fn, 'w') as pdb_file:
            for line in PDB.gen_pdbline(atoms):
                pdb_file.write(line + "\n")
