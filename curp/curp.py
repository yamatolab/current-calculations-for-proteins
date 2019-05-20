"""curp - Main curp programm.

Developed by Yamato's lab, Nagoya University.
Use the argument -h for a list of available arguments.
Input: .cfg file (see http://www.comp-biophys.com/resources/curp-tutorial/tutorials.html
and the tutorial directory for more informations.)
"""

from __future__ import print_function

# standard module
import sys, os
import exception
import time

# curp module
src_dir = os.path.dirname(__file__)
sys.path.insert(0, src_dir)
import clog as logger

# make the process name the CURP.
from setproctitle import setproctitle
setproctitle('curp')

def parse_options():
    """
    Defines how arguments given to curp.py are treated.
    To see optional arguments, type in the terminal python curp.py -h 
    Returns parsed_args, containing all arguments.
    """
    # initialize
    import argparse
    parser = argparse.ArgumentParser(description=
            "Curp program.")

    # definitions
    # parser.add_argument('-a', '--add',
    #         dest='a', type=bool, default=False,
    #         action="store_true",
    #         help='The sample argument.')

    # parser.add_argument('-v', '--vervose',
    #         dest='vervose', type=bool, default=False,
    #         action="store_true",
    #         help='The sample argument.')

    parser.add_argument('-v', '--vervose',
            dest='vervose', default=False,
            action="store_true",
            help='print out informations.')

    parser.add_argument('-s', '--enable-serial',
            dest='use_serial', default=False,
            action="store_true",
            help="calculate in serial, don't calculate in parallel.")

    parser.add_argument('--output-conf-default',
            dest='output_conf_def', default=False,
            action="store_true",
            help=("Output the config parameters in ini format "
                  "with default values."))

    parser.add_argument('--output-conf-formatted',
            dest='output_conf_fmtd', default=False,
            action="store_true",
            help=("Output the config parameters in rest style "
                  "with default values."))

    parser.add_argument(
            'input_', nargs='?', default='run.cfg', #action='append',
            help='specify input filenames.')

    # parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    # parse
    parsed_args = parser.parse_args()
    return parsed_args

"""
Write_ functions below are used to print specific messages and data in the log.
Based mostly on logger.py functions, designed for treating the output which will
appear in the log.
"""


def write_array(array, num_per_line=10):
    """
    print the array (1D tuple, list or numpy array) so that:
    - one line contains ten elements (or less at the end),
    - last character of each element is separated by 5 characters
      from the former element (filled by spaces if needed).
      (Ex: write_array([0.001,0.02,0.1] prints "0.001 0.02  0.1"))
    """

    n = len(array)
    if n/num_per_line != 0:

        for iline_1 in range(n/num_per_line):
            icol_beg, icol_end = 10*iline_1, 10*(iline_1+1)
            line = ' '.join( '{:>5}'.format(col) 
                    for col in array[icol_beg:icol_end] )
            logger.info( line ) #Printed as information.
                                #see clog.py 

        # last line
        if n%num_per_line != 0:
            icol_beg = 10*(iline_1+1)
            line = ' '.join( '{:>5}'.format(col)
                    for col in array[icol_beg:] )
            logger.info( line )

    else:
        # last line
        line = ' '.join( '{:>5}'.format(col) for col in array )
        logger.info( line )


def write_times(label_time_pairs, nstep):
    logger.info_title('The summary of the elasped times')

    msg = '{:>30} : {:12.5f} [s]'

    for i, (label, dt) in enumerate(label_time_pairs):
        if label == 'End':
            # total time
            total_time = sum( time for label, time in label_time_pairs )
            logger.info(msg.format('Total curp time', total_time))
            break

        # elif label == 'Post processing':
            # # average time
            # avg_time = dt / float(nstep)
            # logger.info('    {:>30} : {:12.5f} [s/step]'.format(
                # 'Average run time', avg_time))

        else:
            logger.info(msg.format(label+' time', dt))

def write_success():
    logger.info_title('CURP finished completely.')

    # write the citation to logger
    # citation_fp = os.path.join(src_dir,'..','lib','citation.txt')
    # with open(citation_fp, 'rb') as citation_file:
        # logger.info(citation_file.read())





def init_current(setting, par):
    """
    Input:
      -  setting: Setting class, defined in setting.py
                  contains the settings from .cfg configuration file.
      -  par: SequentialProcessor or ParallelProcessor, cf parallel.py.

    Output:
      -  cal: EnergyFluxCalculator or StressCurrentCalculator class
              defined in current/flux.py or current/current.py.
              made using:
               - interaction table between atoms is made in "if do_table",
                 taking into account the group_pair table (if do_gpair)x.
               - group table made in "if do_group"
               - setting, topology and target_atoms.
      -  writer: MultiFluxWriter or MultiCurrentWriter class,
                 defined in current/writer.py.
                 Used for parallel processes.
      -  topology: TopologyParser class,
                   defined in parser/format/topology.py.
                   Obtained from [input_format] topology_file.
      -  target_atoms: List of Integers,
                       obtained from [curp] target_atoms.
                       Contains the indices of each atom to be considered.
      -  label_time_pairs: List of Tuples (action, time of action),
                           records how many time has elapsed since the
                           beginning of the program and the action.
    """ 

    do_topology   = True
    do_target     = True
    do_group      = True
    do_gpair      = True
    do_table      = True
    do_calculator = True
    do_writer     = True

    label_time_pairs = []

    # confirm whether the specified potential can use the specified format
    # in the input section.

    if do_topology:
        t0 = time.time()

        # get topology
        import parser
        file_format = setting.input.format
        fmt_section = getattr(setting, 'input_' + file_format)
        potential   = setting.curp.potential

        # atom type for stress tensor calculation
        use_atype = True if setting.curp.method == 'momentum-current' else False
        topology = parser.get_tplprm(
                file_format, fmt_section, potential, use_atype)
        """
        topology: TopologyParser class
                   defined in parser/file_format/topology.py
        """
        natom = topology.get_natom()

        # get the list that decompose all potential
        decomp_list = topology.get_decomp_list()

        label_time_pairs += [('Topology', time.time()-t0)]

    if do_target:
        t0 = time.time()

        # parse and get target_atoms
        import table.target as target
        target_atoms = target.parse_target_atoms_line(
                setting.curp.target_atoms, natom)

        logger.info_title('Target atoms information')
        logger.info('target atoms list: ' + ', '
                .join(setting.curp.target_atoms))
        logger.info(' ===> ')
        logger.info('show all atoms explicitly:')
        write_array(target_atoms)
            
        anames = topology.get_atom_info()['names']

        # for atom names
        target_anames = [ str(iatm) + '_' + anames[iatm-1]
                          for iatm in target_atoms ]

        label_time_pairs += [('Target', time.time()-t0)]

    if do_group:
        t0 = time.time()

        # get group atoms by setting.curp.group_method
        from table import group
        
        gname_iatoms_pairs = group.get_group_iatoms_pairs(setting, target_atoms,
                topology.get_residue_info(), topology.get_atom_info() )
        gnames = [ gname for gname, iatoms in gname_iatoms_pairs ]

        label_time_pairs += [('Group', time.time()-t0)]

    if do_gpair:
        t0 = time.time()

        # get group pair table
        # table = {'00001_ALA':['00002_GLY','00003_ASP', ...], ...}
        if setting.curp.group_pair_file[0]:
            from table import group_pair
            gpair_parser = group_pair.GroupPairParser(
                    setting.curp.group_pair_file[0], gnames)
            gpair_table = gpair_parser.get_gpair_table()

        else:
            gpair_table = None

        label_time_pairs += [('Group pair table', time.time()-t0)]

    if do_table:
        t0 = time.time()

        # make table for nonbonded
        nonbonded_table = topology.get_nonbonded_table()

        # make table for target atoms
        if setting.curp.method == 'momentum-current':
            inttable = target.make_interaction_table_current(
                    nonbonded_table, target_atoms, natom)

        elif setting.curp.method in ('energy-flux', 'heat-flux'):
            inttable = target.make_interaction_table_flux(
                    nonbonded_table, target_atoms, natom)

        else:
            inttable = nonbonded_table

        # make table for group pair
        if gpair_table:
            gp = group_pair.GroupPair(gpair_table, gname_iatoms_pairs, natom)
            inttable = gp.get_inttable_with_gpair(inttable)

        # generate interaction table
        interact_table = list(inttable.gen_divided_table())
        memories = inttable.get_table_memories()

        # output interact_table
        logger.info_title('Interaction table informations')
        logger.info('(iatm, jatm_beg, jatm_end)')
        for itab_1, (table, mem) in enumerate(zip(interact_table, memories)):
            msg = '** table {}, memory = {} MB **'
            logger.info( msg.format(itab_1+1, mem/10.0**6) )
            for t in table:
                logger.info(t)

        label_time_pairs += [('Interaction table', time.time()-t0)]

    if do_calculator:
        t0 = time.time()

        # decide and prepare a calculator
        import current
        Calculator = current.get_calculator(setting)
        cal = Calculator()
        cal.prepare(topology=topology, setting=setting,
                target_atoms=target_atoms,
                gname_iatoms_pairs=gname_iatoms_pairs,
                interact_table=interact_table)

        label_time_pairs += [('Calculator setting', time.time()-t0)]

    if do_writer:
        t0 = time.time()

        # decide writer
        if par.is_root():
            from current.writer import get_writer
            writer = get_writer(setting, decomp_list,
                            target_anames, cal.get_groupnames(), gpair_table)
        else:
            writer = None

        label_time_pairs += [('Writer setting', time.time()-t0)]

    return cal, writer, topology, target_atoms, label_time_pairs

def init_dynamics(setting, par):

    do_topology   = True
    do_target     = True
    do_table      = True
    do_calculator = True
    do_writer     = True

    label_time_pairs = []

    # confirm whether the specified potential can use the specified format
    # in the input section.

    if do_topology:
        t0 = time.time()

        # get topology
        import parser
        file_format = setting.input.format
        fmt_section = getattr(setting, 'input_' + file_format)
        potential   = setting.curp.potential

        # atom type for stress tensor calculation
        topology = parser.get_tplprm(
                file_format, fmt_section, potential, use_atomtype=False)
        natom = topology.get_natom()

        label_time_pairs += [('Topology', time.time()-t0)]

    if do_target:
        t0 = time.time()

        # parse and get target_atoms
        import table.target as target
        target_atoms = target.parse_target_atoms_line(
                setting.curp.target_atoms, natom)

        logger.info_title('Target atoms information')
        logger.info('target atoms list: ' + ', '
                .join(setting.curp.target_atoms))
        logger.info(' ===> ')
        logger.info('show all atoms explicitly:')
        write_array(target_atoms)
            
        anames = topology.get_atom_info()['names']

        # for atom names
        target_anames = [ str(iatm) + '_' + anames[iatm-1]
                          for iatm in target_atoms ]

        label_time_pairs += [('Target', time.time()-t0)]

    if do_table:
        t0 = time.time()

        # make table for nonbonded
        nonbonded_table = topology.get_nonbonded_table()
        inttable = nonbonded_table

        # make table for group pair
        # if gpair_table:
            # gp = group_pair.GroupPair(gpair_table, gname_iatoms_pairs, natom)
            # inttable = gp.get_inttable_with_gpair(inttable)

        # generate interaction table
        interact_table = list(inttable.gen_divided_table())
        memories = inttable.get_table_memories()

        # output interact_table
        logger.info_title('Interaction table informations')
        logger.info('(iatm, jatm_beg, jatm_end)')
        for itab_1, (table, mem) in enumerate(zip(interact_table, memories)):
            msg = '** table {}, memory = {} MB **'
            logger.info( msg.format(itab_1+1, mem/10.0**6) )
            for t in table:
                logger.info(t)

        label_time_pairs += [('Interaction table', time.time()-t0)]

    if do_calculator:
        t0 = time.time()

        # get potential function

        # get integrator
        import dynamics
        cal = dynamics.Integrator(topology, setting, interact_table)

        label_time_pairs += [('Calculator setting', time.time()-t0)]

    if do_writer:
        t0 = time.time()

        # decide writer
        # In current version the format of the file to write trajectory can be 
        # only amber format.
        # We are going to allow us select the trajectory format.

        if par.is_root():
            import dynamics
            writer = dynamics.Writer(setting.dynamics)
        else:
            writer = None

        label_time_pairs += [('Writer setting', time.time()-t0)]

    return cal,writer,target_atoms,topology,label_time_pairs

def get_data_iter(setting, topology, target_atoms):
    """
    Input: setting object, topology object, target_atoms object.
    Output: trajectory as an iterable, using "yield" statement.
    Each frame of the trajectory is an element of the iterator.
    """

    # get topology informations
    natom = topology.get_natom()
    masses = topology.get_atom_info()['masses']
    # get flag for beriodic boundary condition
    use_pbc = True if topology.get_pbc_info() else False

    # format
    file_format = setting.input.format
    fmt_section = getattr(setting, 'input_' + file_format)

    # judge wether the potential used in calculation is valid or not.
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
    
    # decide trajectory format and parse trajectory files or restart file
    import parser
    if use_classic:

        data_iter = parser.gen_phase_trajectory(file_format, fmt_section,
                setting.input.first_last_interval, natom, use_pbc,
                setting.curp.remove_trans, setting.curp.remove_rotate,
                target_atoms, masses, logger)

        data_iter.next()
    else:

        data_iter = parser.gen_matrix_ensemple(file_format, fmt_section,
                setting.input.first_last_interval, natom, logger)

    return data_iter


def curp(input_="run.cfg",
         use_serial=False,
         vervose=False,
         output_conf_def=False,
         output_conf_fmtd=False
         ):
    """
    If do_xxx statements are here on readability and debugging purpose
    (steps can be skipped by setting do_xxx as False).

    if do_parallel:
        variable par is created.
        par is a SequentialProcessor() or ParallelProcessor() object (see parallel.py),
        depending on the command line options and if mpi is run or not.
        clog module is configured, setting the log options.

    if do_write:
        Writes informations like date, license or if curp is run in parallel in log.

    if do_setting:
        setting variable is defined as a setting object (see setting.py),
        storing the configuration file informations.
        From there is either do_current or do_dynamics set as True.

    if do_init and do_current:
        init_current() is launched, defining topology or calculator variables.
        par_iter() is launched, defining data_iter, an iterator going through the steps
            of the trajectory.
        par.run() is launched, defining the results_iter as an iterator.

    if do_init and do_dynamics:
        xxxxxxxx

    if do_write:
        results_iter is iterated through, step by step. The calculus starts there.
        Each X step the result is printed in the log.

    if do_summary:
        the job is completed. A summary of the time each process took and such informations
            is printed in the log.
        
    """


    do_parallel = True
    do_title    = True
    do_setting  = True
    do_init     = True
    do_run      = True
    do_write    = True
    do_summary  = True
    if do_parallel:
        import parallel
        if use_serial:
            par = parallel.SequentialProcessor()

        else:
            if parallel.use_mpi:
                par = parallel.ParallelProcessor()
            else:
                par = parallel.SequentialProcessor()

    # configuration of clog module
    if vervose:
        logger.log_level = logger.DEBUG
    else:
        logger.log_level = logger.INFO

    # logger.log_level = logger.WARN

    logger.print_function = par.write

    if output_conf_def:
        import setting as st
        config = st.parse_config()
        setting = st.Setting(config)
        logger.info(setting)
        quit()

    if output_conf_fmtd:
        import setting as st
        config = st.parse_config()
        setting = st.Setting(config)
        logger.info(format(setting, "rst"))
        quit()

    if do_title:
    # Show the Curp's title
    # curp_title = (10*' ' + 60*'=' + 10*' ' + '\n'
            # + 25*' ' + '{:<50}' + 5*' ' + '\n' 
            # + 10*' ' + 60*'-' + 10*' ' + '\n' )
        curp_title = 3*'{:^80}\n'
        logger.info(curp_title.format(60*'=', 'Curp Program Start !! ', 60*'-'))

        # display today's information.
        import datetime
        now = datetime.datetime.today()
        time_fmt = ('TIME STAMP: {year}/{month}/{day:02} '
                + '{hour:02}:{minute:02}:{second:02}' )
        time_stamp = time_fmt.format(year=now.year, month=now.month,
                day=now.day, hour=now.hour, minute=now.minute,
                second=now.second)
        logger.info('{:>80}'.format(time_stamp))

        # display the license content.
        logger.info()
        logger.info()
        license_fp = os.path.join(src_dir, 'LICENSE-short.txt')
        with open(license_fp, 'rb') as license_file:
            for line in license_file:
                logger.info(' '*8 + line.strip())
        logger.info()
        logger.info()

        # display whether the CURP used is serial version or parallel version.
        logger.info_title("Parallel Processing information")
        par_string = 'Use **PARALLEL** caluculation for curp.'
        ser_string = 'Use **SERIAL** caluculation for curp.'
        if use_serial:
            logger.info('{:^80}'.format(ser_string) )
        else:
            if parallel.use_mpi:
                logger.info('{:^80}'.format(par_string) )
            else:
                logger.info('{:^80}'.format(ser_string) )

    if do_setting:
        t0 = time.time()

        # input file
        config_filename = input_

        # # parse and parse setting
        import setting as st
        config = st.parse_config(config_filename)
        setting = st.Setting(config)

        # set the frequency to output log from the setting parameter.
        logger.set_log_frequency(setting.curp.log_frequency)

        label_time_pairs = [('Setting', time.time()-t0)]

    # decide the calculation method
    do_current  = setting.curp.method in ('momentum-current', 'energy-flux', 'heat-flux')
    do_dynamics = setting.curp.method == 'microcanonical'

    if do_init and do_current:
        # log
        logger.info("{:^80}".format("**Current** calculation is performed."))

        # initialize of curp to get the object for calculation
        cal, writer, tpl, target_atoms, label_time_pairs_local = \
                init_current(setting, par)
        label_time_pairs += label_time_pairs_local

        t0 = time.time()
        # write header of data
        if par.is_root(): writer.write_header()

        # get data iterator, ex) trajectory, density matrix, ...
        if par.is_root():
            data_iter = get_data_iter(setting, tpl, target_atoms)
        else:
            data_iter = None
        label_time_pairs += [('Data object parse', time.time()-t0)]

        if do_run:
            t0 = time.time()

            # label_time_pairs += [('Passing to parallel', time.time())]
            ####################################################################
            results_iter = par.run(cal.run, data=data_iter)
            ####################################################################

    elif do_init and do_dynamics:
        logger.info("{:^80}".format("**Dynamics** calculation is performed."))

        # initialize of curp to get the object for calculation
        cal, writer, target_atoms, tpl, label_time_pairs_local \
                = init_dynamics(setting, par)
        label_time_pairs += label_time_pairs_local
        
        t0 = time.time()
        # write header of data
        if par.is_root(): writer.write_header()

        # get data iterator, ex) trajectory
        t0 = time.time()
        if par.is_root():
            data_iter = get_data_iter(setting, tpl, target_atoms)
        else:
            data_iter = None
        label_time_pairs += [('Data object parse', time.time()-t0)]

        if do_run:
            t0 = time.time()

            # label_time_pairs += [('Passing to parallel', time.time())]
            ####################################################################
            results_iter = cal.run(data_iter.next())
            ####################################################################

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
            t1 = time.time()
            cur_steps.append( cur_step )
            logger.debug_cycle('    writing the data at step {} ...'
                    .format(cur_step))
            if par.is_root():
                t_write = time.time()
                writer.write(istep, *results)
                logger.debug_cycle( 'writing at {}: {}'
                        .format(cur_step, time.time()-t_write))
            t2 = time.time()
            dt_writing += t2 - t1

        if par.is_root():
            dt_cal = time.time() - t_cal

            logger.debug_cycle('total/1process : {} / {}'.format(dt_cal,
                dt_cal/float(istep+1)))
                # dt_cal/float(istep+1)/float(par.get_nproc())))

    if do_summary:

        # write the time of writing the data
        label_time_pairs += [('Run', time.time() - t0 - dt_writing )]
        label_time_pairs += [('Writing', dt_writing)]
        # logger.info(msg.format('Writing', dt_writing/float(nstep) ))

        t0 = time.time()

        nstep = istep + 1

        if not par.is_root(): return None

        # if par.is_root():

        if istep == -1:
            msg = "The trajectory data wasn't read"
            raise exception.CurpException(msg)

        msg = '{:>20} time : {:12.5f} [s/snapshot]'

        logger.info_title('Detailed timing in calculator object')
        tag_times_pair_iter = cal.gen_time_info()
        total_time = 0.0
        for tag, times in tag_times_pair_iter:
            each_time = sum(times) / float(len(times))
            logger.info(msg.format(tag, each_time ))
            total_time += each_time
        logger.info('    ' + 49*'=')
        logger.info(msg.format('Total', total_time ))

        # write timing information of parallel run
        logger.info_title('Detailed timing for parallel processing')
        msg = '{:>20} time : {:12.5f} [s/snapshot]'
        time_info = par.get_time_info(nstep)
        for label, t in time_info:
            logger.info(msg.format(label, t ))

        if par.get_other_time():
            msg = '{:>20} time : {:12.5f} [s/all]'
            logger.info(msg.format('misc', par.get_other_time()))
            logger.info(msg.format('Total MPI', par.get_mpi_time()))

        # write the calculated trajectory
        logger.info_title('Steps calculated in trajectory')
        write_array(cur_steps)

        label_time_pairs += [('Post processing', time.time()-t0)]
        label_time_pairs += [('End', 0.0)]

        # write summary
        write_times(label_time_pairs, istep+1)

        # write success
        write_success()

def main():
    # parse command line options
    options = parse_options()
    curp(**vars(options))
    

################################################################################

if __name__ == '__main__':
    main()

    # line = ['1-']

    # print(parse_target_atoms_line(line, 33))
