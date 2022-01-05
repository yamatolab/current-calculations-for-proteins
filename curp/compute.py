"""curp - Main curp program.

Developed by Yamato's lab, Nagoya University.
Use the argument -h for a list of available arguments.
Input: .cfg file
(see https://curp.jp/tutorials/
and the tutorial directory for more informations.)
"""

from __future__ import print_function

# Standard module
import sys
import os
import time

import curp.exception as exception
import curp.clog as logger
import curp.setting as st

# Curp module
SRC_DIR = os.path.dirname(__file__)
sys.path.insert(0, SRC_DIR)

# Write_ functions below are used to print specific messages and data
# in the log.
# Based mostly on logger.py functions, designed for treating the output
# which will appear in the log.


def write_array(array, num_per_line=10):
    """
    Print array (1D tuple, list or numpy array)

    Each line contains num_per_line elements (or less at the end),
    Last character of each element is separated by 5 characters
    from the former element (filled by spaces if needed).

    Parameters
    ----------
    array : array_like
    num_per_line : int, optional
        Number of elements per line

    Examples
    --------
    >>> write_array([0.001, 0.02, 0.1])
    0.001 0.02  0.1
    """

    length = len(array)
    if length//num_per_line != 0:

        for iline_1 in range(length/num_per_line):
            icol_beg, icol_end = 10*iline_1, 10*(iline_1+1)
            line = ' '.join('{:>5}'.format(col)
                            for col in array[icol_beg:icol_end])
            logger.info(line)   # Printed as information, see clog.py

        # Last line
        if length % num_per_line != 0:
            icol_beg = 10 * (length//num_per_line+1)
            line = ' '.join('{:>5}'.format(col)
                            for col in array[icol_beg:])
            logger.info(line)

    else:
        # Last line
        line = ' '.join('{:>5}'.format(col) for col in array)
        logger.info(line)


def write_times(label_time_pairs):
    """Print formatted time with label

    Used to print the time each process took.
    Parameters
    ----------
    label_time_pairs: list of list
        Duration of each process.
        Each element is [label, time].
        label is the process name.
        time is the time the process took.
    """

    logger.info_title('The summary of the elasped times')

    msg = '{:>30} : {:12.5f} [s]'

    for i, (label, dt) in enumerate(label_time_pairs):
        if label == 'End':
            # Total time
            total_time = sum(time for label, time in label_time_pairs)
            logger.info(msg.format('Total curp time', total_time))
            break
        else:
            logger.info(msg.format(label+' time', dt))


def write_success():
    """Print success message."""
    # Write the citation to logger.
    logger.info_title('CURP finished completely.')


def do_topology(setting):
    """Get topology object.

    Parameters
    ----------
    setting : curp.setting.Setting
        Topology file and potential container.

    Returns
    -------
    topology : curp.parser.file_format.topology.TopologyParser
    natom : int
        Number of atom
    """

    import curp.parser as parser
    file_format = setting.input.format
    fmt_section = getattr(setting, 'input_' + file_format)
    potential = setting.curp.potential

    # Atom type depends on stress tensor or flux calculation.
    use_atype = setting.curp.method == 'momentum-current'
    topology = parser.get_tplprm(file_format, fmt_section,
                                 potential, use_atype)
    natom = topology.get_natom()
    return(topology, natom)


def do_target(atom_range, natom):
    """Get target atoms

    Parameters
    ----------
    atom_range : list
        target_atoms field in .cfg file
    topology : curp.parser.file_format.topology.TopologyParser
    natom : int
        Number of atom

    Returns
    -------
    target_atoms : list of int
        List of targeted atoms indexes.
    """

    from curp.table import target
    target_atoms = target.parse_target_atoms_line(atom_range,
                                                  natom)

    logger.info_title('Target atoms information')
    logger.info('target atoms list: ' + ', '.join(atom_range))
    logger.info(' ===> ')
    logger.info('show all atoms explicitly:')
    write_array(target_atoms)

    return target_atoms


def gen_tables(inttable):
    """Generate interaction table.

    Parameters
    ----------
    inttable : table.interact_table.InteractionTable

    Returns
    -------
    interact_table : list
    """

    interact_table = list(inttable.gen_divided_table())
    memories = inttable.get_table_memories()

    # Output interact_table.
    logger.info_title('Interaction table informations')
    logger.info('(iatm, jatm_beg, jatm_end)')
    for itab_1, (table, mem) in enumerate(zip(interact_table, memories)):
        msg = '** table {}, memory = {} MB **'
        logger.info(msg.format(itab_1+1, mem/10.0**6))
        for tab in table:
            logger.info(tab)

    return interact_table


def init_current(setting, par):
    """Initialize current calculator and variables.

    Initialize current, heat flux and energy flux calculations.

    Parameters
    ----------
    setting : curp.setting.Setting
        settings from .cfg configuration file container.
    par : parallel.SequentialProcessor or parallel.ParallelProcessor

    Return
    ------
    cal : curp.current.flux.EnergyFluxCalculator \
          or curp.current.current.StressCurrentCalculator
        Made using
         - interaction table between atoms is made in "if do_table",
           taking into account the group_pair table (if do_gpair).
         - group table made in "if do_group"
         - setting, topology and target_atoms.
    writer : curp.current.writer.MultiFluxWriter or .MultiCurrentWriter
    topology : parser.format.topology.TopologyParser
        Obtained from [input_format] topology_file.
    target_atoms: list of int
        Obtained from [curp] target_atoms.
        Contains the indices of each atom to be considered.
    label_time_pairs: list of tuples (action, time of action),
        Records time taken by each step.
    """

    label_time_pairs = []

    # Get topology object.
    t_0 = time.time()
    topology, natom = do_topology(setting)
    # Get the list that decomposes all potential.
    decomp_list = topology.get_decomp_list()
    label_time_pairs += [('Topology', time.time()-t_0)]

    # Determine target atoms.
    t_0 = time.time()
    target_atoms = do_target(setting.curp.target_atoms, natom)
    # Get target atoms names.
    atom_names = topology.get_atom_info()['names']
    target_anames = [str(i_atom) + '_' + atom_names[i_atom-1]
                     for i_atom in target_atoms]
    label_time_pairs += [('Target', time.time()-t_0)]

    # Get atom groups depending on setting.curp.group_method.
    t_0 = time.time()
    from curp.table import group

    gname_iatoms_pairs = group.get_group_iatoms_pairs(
        setting, target_atoms,
        topology.get_residue_info(),
        topology.get_atom_info())
    gnames = [gname for gname, iatoms in gname_iatoms_pairs]
    label_time_pairs += [('Group', time.time()-t_0)]

    # Get which group interactions will be considered as group pairs.
    t_0 = time.time()
    if setting.curp.group_pair_file[0]:
        from curp.table import group_pair
        gpair_parser = group_pair.GroupPairParser(
            setting.curp.group_pair_file[0], gnames)
        gpair_table = gpair_parser.get_gpair_table()
    else:
        gpair_table = None
    label_time_pairs += [('Group pair table', time.time()-t_0)]

    # Make interaction tables.
    t_0 = time.time()
    # Make table for nonbonded.
    nonbonded_table = topology.get_nonbonded_table()

    # Make table for target atoms.
    from curp.table import target
    if setting.curp.method == 'momentum-current':
        inttable = target.make_interaction_table_current(
            nonbonded_table, target_atoms, natom)
    elif setting.curp.method in ('energy-flux', 'heat-flux'):
        inttable = target.make_interaction_table_flux(
            nonbonded_table, target_atoms, natom)
    else:
        inttable = nonbonded_table

    # Make table for group pair.
    if gpair_table:
        gp = group_pair.GroupPair(gpair_table, gname_iatoms_pairs, natom)
        inttable = gp.get_inttable_with_gpair(inttable)

    interact_table = gen_tables(inttable)

    label_time_pairs += [('Interaction table', time.time()-t_0)]

    # Decide and prepare a calculator.
    t_0 = time.time()

    import curp.current as current
    calculator = current.get_calculator(setting)
    cal = calculator()
    cal.prepare(topology=topology, setting=setting,
                target_atoms=target_atoms,
                gname_iatoms_pairs=gname_iatoms_pairs,
                interact_table=interact_table)
    label_time_pairs += [('calculator setting', time.time()-t_0)]

    # Decide Writer.
    t_0 = time.time()

    if par.is_root():
        from curp.current.writer import get_writer
        writer = get_writer(setting, decomp_list,
                            target_anames, cal.get_groupnames(),
                            gpair_table)
    else:
        writer = None

    label_time_pairs += [('Writer setting', time.time()-t_0)]

    return cal, writer, topology, target_atoms, label_time_pairs


def init_dynamics(setting, par):
    """Initialize microcanonical calculation objects (DEPRECATED).

    Parameters
    ----------
    setting : curp.setting.Setting
        Settings from .cfg configuration file container.
    par : parallel.SequentialProcessor or parallel.ParallelProcessor

    Output:
    cal : curp.current.flux.EnergyFluxCalculat
          or curp.current.current.StressCurrentCalculator
        Made using
         - interaction table between atoms is made in "if do_table",
           taking into account the group_pair table (if do_gpair)x.
         - group table made in "if do_group"
         - setting, topology and target_atoms.
    writer : curp.current.writer.MultiFluxWriter or .MultiCurrentWriter
    topology : parser.format.topology.TopologyParser
        Obtained from [input_format] topology_file.
    target_atoms: list of int
        Obtained from [curp] target_atoms.
        Contains the indices of each atom to be considered.
    label_time_pairs: list of tuples (action, time of action),
        Records time taken by each step.
    """
    label_time_pairs = []

    # Get topology object.
    t_0 = time.time()
    topology, natom = do_topology(setting)

    label_time_pairs += [('Topology', time.time()-t_0)]

    # Determine target atoms.
    t_0 = time.time()
    target_atoms = do_target(setting.curp.target_atoms, natom)
    label_time_pairs += [('Target', time.time()-t_0)]

    # Make table for nonbonded.
    t_0 = time.time()

    inttable = topology.get_nonbonded_table()

    interact_table = gen_tables(inttable)
    label_time_pairs += [('Interaction table', time.time()-t_0)]

    # Defines calculator.
    t_0 = time.time()

    import curp.dynamics as dynamics
    cal = dynamics.Integrator(topology, setting, interact_table)
    label_time_pairs += [('calculator setting', time.time()-t_0)]

    # Decide writer.
    # In current version the format of the file to write trajectory can be
    # only amber format.
    # We are going to allow the selection of trajectory format.
    t_0 = time.time()

    if par.is_root():
        writer = dynamics.Writer(setting.dynamics)
    else:
        writer = None
    label_time_pairs += [('Writer setting', time.time()-t_0)]

    return cal, writer, target_atoms, topology, label_time_pairs


def get_data_iter(setting, topology, target_atoms):
    """Import targeted atoms trajectories

    Parameters
    ----------
    setting : curp.setting.Setting
        Settings from .cfg configuration file container.
    topology : parser.format.topology.TopologyParser
        Topology object obtained from [input_format] topology_file.
    target_atoms: list of int
        Indices of each atom to be considered.
        Obtained from [curp] target_atoms.

    Returns
    -------
    data_iter : generator of tuple
        tuple (cstep_crd, (new_crd, new_vel, pbc_crd))
        cstep_crd : int
            Step of shown coordinates and velocity.
        new_crd : numpy.ma.MaskedArray
            Coordinates of each atom at step cstep_crd.
        new_vel : numpy.ma.MaskedArray
            Velocity of each atom at step cstep_crd.
        pbc_crd : None or tuple of float
            Periodic boundary condition coordinates.
    """

    # Get topology informations.
    natom = topology.get_natom()
    masses = topology.get_atom_info()['masses']
    # Get flag for periodic boundary condition.
    use_pbc = True if topology.get_pbc_info() else False

    # Format
    file_format = setting.input.format
    fmt_section = getattr(setting, 'input_' + file_format)

    # Judge wether the potential used in calculation is valid or not.
    method = setting.curp.method
    if method == 'momentum-current':
        use_classic = True

    elif method == 'energy-flux':
        use_classic = True

    elif method == 'heat-flux':
        use_classic = True

    elif method == 'electron-transfer':
        use_classic = False

    elif method == 'microcanonical':
        use_classic = True

    else:
        pass

    # Decide trajectory format and parse trajectory files or restart file.
    import curp.parser as parser

    if use_classic:
        data_iter = parser.gen_phase_trajectory(
            file_format, fmt_section, setting.input.first_last_interval,
            natom, use_pbc, setting.curp.remove_trans,
            setting.curp.remove_rotate, target_atoms, masses, logger)

        data_iter.next()
    else:
        data_iter = parser.gen_matrix_ensemble(
            file_format, fmt_section, setting.input.first_last_interval,
            natom, logger)

    return data_iter


def curp(input_="run.cfg", use_serial=False, vervose=False,
         output_conf_def=False, output_conf_fmtd=False, **kwds):
    """Compute stress tensor or inter-residue energy or heat flow.

    The computation steps are as follow:
    if do_parallel:
        variable par is created.
        par is a SequentialProcessor() or ParallelProcessor() object
        (see parallel.py),
        depending on the command line options and if mpi is run or not.
        clog module is configured, setting the log options.
    if do_write:
        Writes informations like date, license or if curp is run in
        parallel in log.
    if do_setting:
        setting variable is defined as a setting object (see setting.py),
        storing the configuration file informations.
        From there is either do_current or do_dynamics set as True.
    if do_init:
        init_current() is launched, defining topology and calculator
        variables.
        par_iter() is launched, defining data_iter, a generator going
        through the steps of the trajectory.
        par.run() is launched, defining the results_iter as a generator.
    if do_write:
        results_iter is iterated through, step by step. The calculus
        starts there.
        Each X step the result is printed in the log.
    if do_summary:
        The job is completed. A summary of the time each process took
        and such informations is printed in the log.

    Parameters
    ----------
    input : list of strings
        Configuration files names.
    use_serial : bool, optional
        Whether computations are run in parallel or not (is serial).
    vervose : bool, optional
        If True prints debug level messages in the log.
    output_conf_def :
        If True outputs the config parameters in ini format with default
        values.
    output_conf_fmtd :
        If True outputs the config parameters in rest style with default
        values.
    """

    do_parallel = True
    do_title = True
    do_setting = True
    do_init = True
    do_run = True
    do_write = True
    do_summary = True
    if do_parallel:
        import curp.parallel as parallel
        if use_serial:
            par = parallel.SequentialProcessor()
        elif parallel.use_mpi:
            par = parallel.ParallelProcessor()
        else:
            par = parallel.SequentialProcessor()

    # Configuration of clog module.
    if vervose:
        logger.log_level = logger.DEBUG
    else:
        logger.log_level = logger.INFO

    logger.print_function = par.write

    if output_conf_def:
        config = st.parse_config()
        setting = st.Setting(config)
        logger.info(setting)
        quit()

    if output_conf_fmtd:
        config = st.parse_config()
        setting = st.Setting(config)
        logger.info(format(setting, "rst"))
        quit()

    # Show Curp's title.
    if do_title:
        curp_title = 3*'{:^80}\n'
        logger.info(curp_title.format(60*'=',
                                      'Curp Program Start !! ',
                                      60*'-'))
        # Display today's information.
        import datetime
        now = datetime.datetime.today()
        time_fmt = ('TIME STAMP: {year}/{month}/{day:02} '
                    + '{hour:02}:{minute:02}:{second:02}')
        time_stamp = time_fmt.format(year=now.year, month=now.month,
                                     day=now.day, hour=now.hour,
                                     minute=now.minute, second=now.second)
        logger.info('{:>80}'.format(time_stamp))

        # Display the license content.
        logger.info()
        logger.info()
        license_fp = os.path.join(SRC_DIR, 'LICENSE-short.txt')
        with open(license_fp, 'rb') as license_file:
            for line in license_file:
                logger.info(' '*8 + line.strip())
        logger.info()
        logger.info()

        # Display whether the CURP used is serial version or parallel version.
        logger.info_title("Parallel Processing information")
        par_string = 'Use **PARALLEL** caluculation for curp.'
        ser_string = 'Use **SERIAL** caluculation for curp.'
        if use_serial:
            logger.info('{:^80}'.format(ser_string))
        else:
            if parallel.use_mpi:
                logger.info('{:^80}'.format(par_string))
            else:
                logger.info('{:^80}'.format(ser_string))

    if do_setting:
        t_0 = time.time()

        # Input file
        config_filename = input_

        # Parse and parse setting
        config = st.parse_config(config_filename)
        setting = st.Setting(config)

        # Set the frequency to output log from the setting parameter.
        logger.set_log_frequency(setting.curp.log_frequency)

        label_time_pairs = [('Setting', time.time()-t_0)]

    # Decide the calculation method
    do_current = setting.curp.method in ('momentum-current',
                                         'energy-flux',
                                         'heat-flux')
    do_dynamics = setting.curp.method == 'microcanonical'

    if do_init and do_current:
        # log
        logger.info("{:^80}".format("**Current** calculation is performed."))

        # Initialize of curp to get the object for calculation.
        cal, writer, tpl, target_atoms, label_time_pairs_local = \
            init_current(setting, par)
        label_time_pairs += label_time_pairs_local

        t_0 = time.time()
        # Write header of data.
        if par.is_root():
            writer.write_header()
            # Get data iterator. Ex) trajectory, density matrix, ...
            data_iter = get_data_iter(setting, tpl, target_atoms)
        else:
            data_iter = None
        label_time_pairs += [('Data object parse', time.time()-t_0)]

        if do_run:
            t_0 = time.time()

            ############################################################
            results_iter = par.run(cal.run, data=data_iter)
            ############################################################

    elif do_init and do_dynamics:
        logger.info("{:^80}".format("**Dynamics** calculation is performed."))

        # Initialize of curp to get the object for calculation.
        cal, writer, target_atoms, tpl, label_time_pairs_local \
            = init_dynamics(setting, par)
        label_time_pairs += label_time_pairs_local

        t_0 = time.time()
        if par.is_root():
            # Write header of data.
            writer.write_header()
            # Get data iterator. Ex: trajectory.
            data_iter = get_data_iter(setting, tpl, target_atoms)
        else:
            data_iter = None
        label_time_pairs += [('Data object parse', time.time()-t_0)]

        if do_run:
            t_0 = time.time()

            ############################################################
            results_iter = cal.run(data_iter.next())
            ############################################################

    else:
        logger.info('ERROR: The method given in the CURP is not available.')
        exit()

    if do_write:

        logger.info_title("The each step's calculated informations")
        cur_steps = []
        istep = -1

        dt_writing = 0.0
        if par.is_root():
            t_cal = time.time()
        for istep, (cur_step, results) in enumerate(results_iter):
            t_1 = time.time()
            cur_steps.append(cur_step)
            logger.debug_cycle('    writing the data at step {} ...'
                               .format(cur_step))
            if par.is_root():
                t_write = time.time()
                writer.write(istep, *results)
                logger.debug_cycle('writing at {}: {}'
                                   .format(cur_step, time.time()-t_write))
            t_2 = time.time()
            dt_writing += t_2 - t_1

        if par.is_root():
            dt_cal = time.time() - t_cal

            if istep > -1:
                logger.debug_cycle('total/1process : {} / {}'
                                   .format(dt_cal, dt_cal/float(istep+1)))

    if do_summary:

        # Write the time of writing the data.
        label_time_pairs += [('Run', time.time() - t_0 - dt_writing)]
        label_time_pairs += [('Writing', dt_writing)]
        t_0 = time.time()

        nstep = istep + 1

        if not par.is_root():
            return None

        if istep == -1:
            msg = "The trajectory data wasn't read"
            raise exception.CurpException(msg)

        msg = '{:>20} time : {:12.5f} [s/snapshot]'

        logger.info_title('Detailed timing in calculator object')
        tag_times_pair_iter = cal.gen_time_info()
        total_time = 0.0
        for tag, times in tag_times_pair_iter:
            each_time = sum(times) / float(len(times))
            logger.info(msg.format(tag, each_time))
            total_time += each_time
        logger.info('    ' + 49*'=')
        logger.info(msg.format('Total', total_time))

        # Write timing information of parallel run.
        logger.info_title('Detailed timing for parallel processing')
        msg = '{:>20} time : {:12.5f} [s/snapshot]'
        time_info = par.get_time_info(nstep)
        for label, t in time_info:
            logger.info(msg.format(label, t))

        if par.get_other_time():
            msg = '{:>20} time : {:12.5f} [s/all]'
            logger.info(msg.format('misc', par.get_other_time()))
            logger.info(msg.format('Total MPI', par.get_mpi_time()))

        # Write the calculated trajectory.
        logger.info_title('Steps calculated in trajectory')
        write_array(cur_steps)

        label_time_pairs += [('Post processing', time.time()-t_0)]
        label_time_pairs += [('End', 0.0)]

        # Write the time each process took and success.
        write_times(label_time_pairs)
        write_success()


########################################################################


if __name__ == '__main__':
    from console import arg_compute, exec_command

    parser = arg_compute()
    exec_command(parser)
