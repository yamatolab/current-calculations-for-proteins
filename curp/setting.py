"""setting - Defines the Setting class

Setting class contains all sections and variables from .cfg
configuration file.
Defines all the .cfg sections as class variables (Section object)
"""

from curp.setting_base import *


class Setting(SettingBase):
    """Settings for curp computations

    Parameters
    ----------
    config : configparser.SafeConfigParser
    check : Bool
        Whether entries are checked or not.

    Examples
    --------
    To create a default Setting object:
    >>> import curp.setting as st
    >>> config = st.parse_config("")
    >>> setting = st.Setting(config)

    To create a Setting object from a configuration file:
    >>> import curp.setting as st
    >>> config = st.parse_config(path_to_cfg)
    >>> setting = st.Setting(config)

    To see the sections:
    >>> dir(setting)

    To print a usable .cfg file
    >>> print(setting)

    To print a usable .rst summary
    >>> print(format(setting, "rst"))

    setting.session_name to access only a section
    Ex: input section
    >>> setting.input

    To show only the values of the input section
    >>> print(setting.input)

    To modify a variable
    Ex: To change curp method to energy-flux
    >>> setting.curp.method = "energy-flux"
    """

    input = Section(

        description = "",

        format = Choice(default='amber', require=True,
            desc='The format of the various files generated by other grograms.',
            values=['presto', 'amber'], value_type=String),

        first_last_interval = List(default=(0, -1, 1), require=False,
             value_type=Int,
             desc=('The first, last and interval step to read coordinates and '
                 'velocities trajectory.'),
             ),

        use_simtime = Bool(default=True, require=False,
            desc='use the simulation time information in trajectory.'),

        group_file = File(default='group.cfg',require=False, allow_glob=False,
            desc=('The group file path to define group.'
                'if curp.group_method == file then this definition is used.')
        ),

        # group_excluded = List(default=[], require=False, value_type=String,
        #     desc='The excluded group list specified in section name'),

    )

    input_amber = Section(

        target = Choice(default='trajectory', require=False,
            desc='The target input file to be used.',
            values=['trajectory', 'restart'], value_type=String),

        topology_file = File(default='', require=True, allow_glob=False,
            desc='Topology file path.'),

        coordinate_file = File(default='', require=False, allow_glob=True,
            desc='Coordinate trajectory file path.'),

        velocity_file = File(default='', require=False, allow_glob=True,
            desc = 'Velocity trajectory file path.'),

        restart_file = File(default='', require=False, allow_glob=False,
            desc = 'Restart file path.'),

        coordinate_format = Choice(default='ascii', require=False,
            desc='The format of coordinate file.',
            values=['ascii', 'netcdf'],
            value_type=String),

        velocity_format = Choice(default='ascii', require=False,
            desc='The format of velocity file.',
            values=['ascii', 'netcdf'],
            value_type=String),

        restart_format = Choice(default='restart', require=False,
            desc='The format of restart file.',
            values=['restart'],
            value_type=String),

        dump_parameters = Bool(default=False, require=False,
            desc='Dump the parsed Amber force field parameter set.'),

    )

    # input_presto = Section(

        # target = Choice(default='trajectory', require=False,
            # desc='The target input file to be used.',
            # values=['trajectory', 'restart'], value_type=String),

        # topology_file = File(default='', require=True, allow_glob=False,
            # desc='Topology file name.'),

        # coordinate_file = File(default='', require=False, allow_glob=True,
            # desc='Coordinate trajectory file name.'),

        # velocity_file = File(default='', require=False, allow_glob=True,
            # desc = 'Velocity trajectory file name.'),

        # restart_file = File(default='', require=False, allow_glob=False,
            # desc = 'Restart trajectory file name.'),

        # coordinate_format = Choice(default='ascii', require=False,
            # desc='The format of coordinate file.',
            # values=['ascii', 'bin'],
            # value_type=String),

        # velocity_format = Choice(default='ascii', require=False,
            # desc='The format of velocity file.',
            # values=['ascii', 'bin'],
            # value_type=String),

        # restart_format = Choice(default='ascii', require=True,
            # desc='The format of restart file.',
            # values=['ascii', 'bin'],
            # value_type=String),

    # )

    volume = Section(

        method = Choice(default='voronoi', require=False,
            desc = 'Algorithm to calculate the atomic volumes.',
            values = ['none', 'vdw', 'voronoi', 'outer'], #, 'smve'],
            value_type=String),

        # smve_rmax = Float(default=3.0, require=False,
        #     desc = ('The cutoff length for the calculation of radial '
        #     'distribution function.')),

        # smve_dr = Float(default=0.01, require=False,
        #     desc = 'The width of bin for radial distribution function.')

        # smve_interval = Int(default=1, require=False,
        #     desc = ('The step interval in reading the trajectory for radial '
        #     'distribution function.')),

        # smve_increment = Int(default=5, require=False,
        #     desc = 'Search the 1st wells every "increment" step'),

        atomic_trajectory_file = File(default='', require=False,
            allow_glob=False,
            desc='Atomic volumes trajectory file path for outer method.'),

        group_trajectory_file = File(default='', require=False,
            allow_glob=False,
            desc='Group volumes trajectory file path for outer method.'),

        voronoi_cutoff = Float(default=6.0, require=False,
            desc = ('The cutoff length that the voronoi calculation finds out'
                    'neighbour candidate particles.') ),

        voronoi_no_hydrogen = Bool(default=False, require=False,
            desc=('Flag to determine whether include hydrogen atoms'
                  'for the voronoi calculation.') ),

        voronoi_solvation = Choice(default='none', require=False,
            desc=('The kind of solvation system to sink the target system '
                  'in vacuum for the voronoi method. '
                  'The density value of the water under NPT ensemble is '
                  '0.99651 [g/cm^3] at 27 [Kelvin]'),
            values=['none', 'RANDOM20', ], #'NPT20', 'NPT25' 'NPT30'],
            value_type=String),

        voronoi_probe_length = Float(default=2.4, require=False,
            desc = ('The probe length of the solvation for the voronoi method.'
                'The water molecules within the probe length '
                'from the system are removed.')
            ),

        voronoi_output_solvation_file = File(default='', require=False,
            allow_glob=False,
            desc=('The file path to write out the solvation pdb data in the '
                  'case of voronoi_solvation == "none".'
                  'If the file path is not given, writing out is not performed.')            ),

        output_volume_file = File(default='', require=False,
            allow_glob=False,
            desc=('The file path to write out the atomic volumes trajectory. '
                  'If this value is not given, writing out is not performed.'
                  'The file written by this option can be used in the options,'
                  'atomic_trajectory_file.'),
            ),

        output_gvolume_file = File(default='', require=False,
            allow_glob=False,
            desc=('The file path to write out the group volumes trajectory. '
                  'If this value is not given, writing out is not performed.'
                  'The file written by this option can be used in the options'
                  'group_trajectory_file.'),
            ),

    )

    curp = Section(

        potential = Choice(default='amberbase', require=True,
            desc='The potential function to calculate the pairwise forces.',
            values=['amberbase', 'amber94', 'amber96', 'amber99', 'amber99SB',
                'amber03', 'amber12SB'], # 'amber-polar'],
            value_type=String),

        target_atoms = List(default=['1-'], require=False, value_type=String,
            desc='The atom list calculated.'),

        method = Choice(default='momentum-current', require=True,
            value_type=String,
            desc=('The method of calculation.'
                '"momentum-current" calculates the stress tensor for systems. '
                '"energy-flux" calculates the energy flow for systems. '
                '"kinetic-flux" calculates the kinetic energy flow for systems. '
                '"dynamics is mainly used to verify the validity of '
                ' the CURP program numerically, '
                'so its implementation is very simple'),
            values=['energy-flux', 'momentum-current',
                    'microcanonical', 'heat-flux', 'kinetic-flux']),
                    #'energy-current', 'stress-flux']),

        group_method = Choice(default='none', require=False,
            value_type=String,
            desc=('The method to construct the group.'
                '"united" means that hydrogen atoms are included in hevy atoms '
                'covalent to them.'
                '"residue" means that the groups are calculated by residue '
                'level. If "file" is specified, the groups definition is given'
                'by the group file in input section.'),
            values=['united', 'residue', 'file', 'none'] ),

        flux_grain = Choice(default='group', require=False,
            value_type=String,
            desc=('The grain to calculate the flux.'
                '"atom", "group" and "both" values mean that the flux '
                'for atom parirs, group pairs and both of them will be '
                'calculated, respectively.'),
            values=['atom', 'group', 'both'] ),

        decomp_group_current = Bool(default=False, require=False,
            desc=('Flag whether decompose the group current into inside and '
                  'outside contributions of group region.'
                  'This option is used for calculating momentum current') ),

        group_pair_file = File(default='',require=False, allow_glob=False,
            desc=("Path to file to define group pair. "
                  "If you didn't given, all of pairs within the targets "
                  "will be calculated.")),

        coulomb_method = Choice(default='cutoff', require=False,
            desc='The method to calculate coulomb interaction.',
            values=['cutoff',], # 'multipole'],
            value_type=String),

        coulomb_cutoff_method = Choice(default='atom', require=False,
            desc='The method to cut off the coulomb interaction.',
            values=['atom'], # 'residue-center','residue-while'],
            value_type=String),

        coulomb_cutoff_length = Float(default=99.9, require=False,
            desc='The cutoff length for the coulomb interaction.'),

        # vdw_method = Choice(default='cutoff', require=False,
        #     desc='The method to calculate van der Waals interaction.',
        #     values=['cutoff', 'multipole'], value_type=String)

        vdw_cutoff_method = Choice(default='atom', require=False,
            desc='The method to cut off the van der Waals interaction.',
            values=['atom'], # 'residue-center','residue-while'],
            value_type=String),

        vdw_cutoff_length = Float(default=99.9, require=False,
            desc='The cutoff length for the van der Waals interaction.'),

        # postproccessing = Choice(default='output', require=False,
        #     desc='The define what the post-proccessing part use.',
        #     values=['output'], # 'residue-center','residue-while'],
        #     value_type=String),

        enable_inverse_pair = Bool(default=False, require=False,
                desc=('Calculate and write out inverse pairs j <- i for flux '
                    'adding normal group pairs: i <- j. '
                    'This option is used in the case '
                    'for calculating energy flux.')),

        remove_trans = Bool(default=True, require=False,
            desc=('Remove the coordinate and velocity of translation '
                  'for the target atoms')),

        remove_rotate = Bool(default=True, require=False,
            desc=('Remove the coordinate and velocity of rotation '
                  'for the target atoms')),

        log_frequency = Int(default=10, require=False,
            desc='Log informations will be written out every given steps.'),

        # restraint
    )

    dynamics = Section(

        integrator = Choice(default='vverlet', require=True,
            desc='The integrator to want to use with the dynamics.',
            values=['vverlet', 'leapfrog'],
            value_type=String),

        dt = Float(default=0.001, require=True,
            desc='Time step to advance snapshots to next step, in ps unit.'),

        num_steps = Int(default=1, require=True,
            desc='The number of integration steps.'),

        crds_file = File(default='', require=False, allow_glob=False,
            desc=("The file path to write out the coordinates trajectory. "
                  "If empty, then don't write.")),

        vels_file = File(default='', require=False, allow_glob=False,
            desc=("The file path to write out the velocities trajectory. "
                  "If empty, then don't write.")),

        crds_frequency = Int(default=1, require=False,
            desc='The interval step to write coordinate trajectory.'),

        vels_frequency = Int(default=1, require=False,
            desc='The interval step to write velocity trajectory.'),

        trj_format = Choice(default='ascii', require=False,
            desc='The format of coordinates and velocities trajectory file.',
            values=['ascii', 'netcdf'],
            value_type=String),
        )

    output = Section(

        # current or flux data
        filename = File(default='current.dat', require=False, allow_glob=False,
            desc='The file name to output the current or flux information.'),

        format = Choice(default='ascii', require=False,
            desc='The format of flux data.',
            values=['ascii', 'netcdf'],
            value_type=String),

        decomp = Bool(default=False, require=False,
            desc=('Flag whether decompose and output the total current/flux,'
                  'autocorrelation function and heat conductivity'
                  'to bonded, coulomb, and van der Waals interaction.') ),

        frequency = Int(default=0, require=False,
            desc=('The frequency to create new file to write '
                'the additional file.')),

        compress = Bool(default=False, require=False,
            desc=('Flag whether compress with gnu zip, '
                'then the extension of the file name became ".gz".') ),

        # basic data
        output_energy = Bool(default=False, require=False,
            desc='Flag whether output the energy information or not.'),

        energy_file = File(default='energy.dat', require=False,
            desc='The file path to output the energy information.'),

        energy_decomp = Bool(default=False, require=False,
            desc=('Flag whether decompose the total energy '
                'to bonded, coulomb, and van der Waals interaction.') ),

        energy_freqency = Int(default=0, require=False,
            desc=('The frequency to write the energy information.')),

        energy_compress = Bool(default=True, require=False,
            desc=('Flag whether compress with gnu zip for energy_file, '
                'then the extension of the file name became ".gz".')),

        # output_force = Bool(default=False, require=False,
            # desc='Flag whether output the force information or not.'),

        # force_file = File(default='force.dat', require=False,
            # desc='The file name to output the force information.'),

        # force_decomp = Bool(default=False, require=False,
            # desc=('Flag whether decompose the total force'
                # 'to bonded, coulomb, and van der Waals interaction.') ),

        # force_freqency = Int(default=0, require=False,
            # desc='The frequency to create new file to write other file.'),

        # force_compress = Bool(default=True, require=False,
            # desc=('Flag whether compress by gnu zip, '
                # 'then the extension of the file name became ".gz".') ),

        # output_tbforce = Bool(default=False, require=False,
            # desc=('Flag whether output the tbforce information or not.'
                # 'note that the file size became too large.')),

        # tbforce_file = File(default='tbforce.dat', require=False,
            # desc='The file name to output the tbforce information.'),

        # tbforce_decomp = Bool(default=False, require=False,
            # desc=('Flag whether decompose the total tbforce'
                # 'to bonded, coulomb, and van der Waals interaction.') ),

        # tbforce_freqency = Int(default=0, require=False,
            # desc='The frequency to create new file to write other file.'),

        # tbforce_compress = Bool(default=True, require=False,
            # desc=('Flag whether compress by gnu zip, '
                # 'then the extension of the file name became ".gz".') )
    )


def parse_config(config_filename=None):
    try:
        import configparser as cp
    except ImportError:
        import ConfigParser as cp
    except:
        raise

    try:
        from cStrinIO import StringIo
    except ImportError:
        from io import StringIO
    except:
        raise

    outfile = StringIO()

    if config_filename:
        with open(config_filename, 'r') as config_file:
            for line in config_file:
                body_ = line.split('#')[0]
                body  = body_.split(';')[0]
                outfile.write(body)

        outfile.seek(0)
    
    import curp.clog as logger
    logger.debug("Type outfile: ", type(outfile))
    logger.debug(outfile)
    
    config = cp.SafeConfigParser()

    config.readfp(outfile)

    outfile.close()

    # config = cp.SafeConfigParser()
    # config.read(config_filename)

    return config

if __name__ == '__main__':
    config = parse_config('run.cfg')
    config = parse_config()
    setting = Setting(config, check=False)
    # setting = Setting(config)
    # print(setting)

    print(format(setting, "rst"))
