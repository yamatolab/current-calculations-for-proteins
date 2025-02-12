=====================
Heat current analysis
=====================

**Local heat current**

The thermal conductivity, :math:`\lambda`, is expressed in terms of the heat current autocorrelation function (HCACF) as:

.. math::
   \lambda = \frac{1}{3Vk_BT^2} \int_0^\infty \left\langle \mathbf{h}(0) \cdot \mathbf{h}(t) \right\rangle dt.

The heat current, :math:`\mathbf{h}`, is written as a summation of interatomic heat current, :math:`\mathbf{h}_{ij}`, as:

.. math::
   \mathbf{h} = \sum_i \sum_{j>i}(\mathbf{r}_i - \mathbf{r}_j)J_{i \leftarrow j} = \sum_i \sum_{j>i} \mathbf{h}_{ij}



The heat current between atom groups, A and B, can be calculated by the summation of local heat current between the pair:

.. math::

   \mathbf{h}_{AB} =
      \sum_{i \in A} \sum_{j \in B} \mathbf{h}_{ij}.

If the atom groups A and B represent residues A and B, respectively, the contribution of this local heat current between the residue pair to the overall thermal conductivity is expressed as:

.. math::

   \Lambda_{AB} = \int_0^\infty  {\left\langle
      {\mathbf{h}_{AB}}(t) \cdot {\mathbf{h}_{AB}}(0)
   \right\rangle dt}.

This calculation is performed by the numerical integration of the auto-correlation function of :math:`\mathbf{h}` using trapezoidal rule.  The CURP program performs these calculations for the proteins and visualizes the energy conductivities using an energy exchange network (EEN). This allows the user to grasp the pattern of the inter-residue interactions at a glance. In principle, it is possible to calculate :math:`\Lambda` based on a single NVE trajectory. However, it is highly recommended to calulate thermal averages of :math:`\Lambda` based on multiple NVE trajectories that start from different points on the energy landscape of the target protein. If the number of the NVE trajectories are sufficiently large, the thermally averaged :math:`\Lambda` are fair representations of the :math:`\Lambda` of the thermally fluctuating protein system.  


CURP analyzes the inter-residue interactions through the following steps:

.. contents::
   :local:
   :depth: 1

Preparing trajectories
=======================

Solvent stripping 
------------------

When we are focusing on residue-residue interactions in a polypeptide chain, solvent coordinates and velocities are not needed. It is possible to save disc space by running the follwoing scripts and disregarding the solvent degrees of freedom:

.. code-block:: bash

   # Solvent spripping from the coordinates trajetory
   curp conv-trj -crd \
       -p system.prmtop   -pf amber \
       -i nev.crd.nc      -if netcdf --irange 1 -1 1 \
       -o stripped.crd.nc -of netcdf --orange 1 -1 1 \
       mask -m mask.pdb

   # Solvent stripping from the velocities trajectory
   curp conv-trj -vel \
       -p system.prmtop   -pf amber \
       -i nev.vel.nc      -if netcdf  --irange 1 -1 1 \
       -o stripped.vel.nc -of netcdf  --orange 1 -1 1 \
       mask -m mask.pdb

Suppose that the time step of the NVE simulation was 2 fs.
In the first (second) example above, the `mask` subcommand removes the solvent degrees of freedom from the coordinates (velocities) trajectory file, `nev.crd.nc` (`nev.vel.nc`), in which coordinates are saved every 10 fs = 5 frames (2 fs = 1 frame), and a solvent-stripped trajectory is stored in `stripped.crd.nc` (`stripped.nev.nc`). The filename of the parameter/topology file is specified after the `-p` option.

**mask**
   This subcommand extracts selected atoms, which are recorded in the PDB file specified after the `-m` option. Make sure that the PDB file contains only the selected atoms, and each atom id coincides with the sequential number of the corresponding atom in the trajectory file. The other atoms will be excluded from the system.

`conv-trj` command options

`-crd`
   Setting coordinates trajectory files, whose type is assumed to be the coordinate type by default.
   Note that this option is mutually exclusive with the `-vel` option.

`-vel`
   Setting velocities trajectories files.
   Note that this option is mutually exclusive with the `-crd` option.

`-p` 
   Specifies the parameter and topology file.

`-pf`
   Specifies the format of the parameter and topology file.

`-i`
   This option specifies the input trajectory file. By repeating this option,
   multiple trajectory files are read in the order you provided.
   
`-if`
   Specifies the format of the input trajectory file.

`\--irange`
   Specifies <first_frame>, <last_frame> and <frame_interval> for the input trajectory file.
   `-1` for <last_frame> represents the last frame of the trajectory.
   For example, if "1 -1 5" is provided, the conv-trj command reads
   1st, 6th, 11th, ... , up to the last frames out from the trajectory.

`-o`
   Specifies the trajectory file for output. You are not allowed to
   specify this option multiple times.

`-of`
   Specifies the format of the output trajectory file.

`\--orange`
   Specifies <first_frame>, <last_frame> and <frame_interval> for the output trajectory file.
   `-1` for <last_frame> represents the last frame of the trajectory.

Note that the parameter and topology file is needed to be modified 
according to the solvent splitting for the subsequent MD simulations of the
new system.

Adjusting the time points for the coordinates and velocities trajectory
------------------------------------------------------------------------

In the Amber restart and trajectory files, the time frames of atomic velocities are shifted by :math:`-\Delta t/2` from those of atomic coordinates.
In the CURP program, however, the time points of the atomic coordinates must coincide with those of the atomic velocities.
Therefore, the atomic velocities in the Amber trajectory file must be 
modified before the energy flow calculations. As explained below, 
the `conv-trj` program of the CURP package estimates the atomic velocities at
the time point of :math:`t` from those at the time points of
:math:`t - \Delta t/2`, and :math:`t + \Delta t/2`. Accordingly, the
velocity frames should be recorded every step to the trajectory file. 

.. note:: Suppose that you started your MD simulation from time :math:`t_0`. Time points of the Amber coordinate and the velocity frames are, then, (:math:`t_0 + \Delta t, t_0 + 2\Delta t, t_0 + 3\Delta t, \cdots`), (:math:`t_0 + \Delta t/2, t_0 + 3\Delta t/2, t_0 + 5\Delta t/2, \cdots`), respectively. In addition, the restart file contains the coordinate and the velocity frames at the time points of :math:`t_0` and :math:`t_0 - \Delta t/2`, respectively.

To modify the time points of the Amber velocities trajectory, the following 
script is available:

.. code-block:: bash

   # adjust the velocity time
   curp conv-trj -vel \
       -p stripped.prmtop -pf amber \
       -i stripped.vel.nc -if netcdf --irange 1 -1 1 \
       -o adjusted.vel.nc -of netcdf --orange 5 -1 5 \
       adjust-vel

Suppose the the NVE simulation was performed with the time step of 2 fs.
In the above example, the 1st, 2nd, ... , up to the last frames are read from 
the `stripped.vel.nc` file obtained in the `Solvent stripping`_ section.
For each of the frame-pairs, (4th, 5th), (9th, 10th), :math:`\cdots` , a new frame is generated at the midtime point of the frame-pair. Consequently, the time points of the the original 5th, 10th, :math:`\cdots`, frames are shifted by :math:`- \Delta t/2` altogether and output to the `adjusted.vel.nc` file. Note that this process is needed for the velocities trajectory obtained by the leap frog algorithm, while not needed for that obtained by the velocity Verlet algorithm.

**adjust-vel**
   This subcommand shifts the time points of the velocities trajectories as described above.  

In this tutorial, we provided the sample trajectories in which the 
time points of the coordinates and velocities trajectories were adjusted.

Avoiding missing frame problem while concatenating multiple trajectory files
----------------------------------------------------------------------------

A special care is needed when you use the leap frog integrator, which is usually employed in the AMBER program, and split your trajectory into multipile files. Suppose that the time period of the `i`-th trajectory segment is :math:`[T_{i-1} + \Delta t, T_{i}]`, and the atomic coordinates and the velocities in this segment are recorded at the time points of :math:`(T_{i-1} + \Delta t, T_{i-1} + 2\Delta t, \cdots ,T_{i} - \Delta t, T_{i})` and :math:`(T_{i-1} + \Delta t/2, T_{i-1} + 3\Delta t/2, \cdots, T_{i} - 3\Delta t/2, T_{i} - \Delta t/2)`, respectively. If you need to consider the velocities at :math:`T_{i}` for the further calculations of energy flow, you need the velocity trajectory files of both `i`-th and `(i+1)`-st segments, because the velocities at :math:`T_{i}` are estimated from those at :math:`T_{i} - \Delta t/2` and :math:`T_{i} + \Delta t/2`. Note that the velocity trajectory file of the `i`-th segment can be replaced with the restart file generated at the end of the `i`-th segment.

**Example**
    When you perform a MD simulation for 10 ps with the time step of :math:`\Delta t` = 2 fs, and save the atomic coordinates every 10 fs, the time points of the atomic coordinates are 10 fs, 20 fs, :math:`\cdots`, 9990 fs, and 10000 fs. On the other hand, you need to save the velocities every step because you need to adjust the time points of the velocities to those of the atomic coordinates. Suppose that you divide the velocities trajectory into halves, and save the first (second) half to the trajectory file named `nve1.vel.nc` (`nve2.vel.nc`). Then the time points of the velocities in `nve1.vel.nc` (`nve2.vel.nc`) are 1 fs, 3 fs, :math:`\cdots`, 4999 fs (5001 fs, 5003 fs, :math:`\cdots`, 9999 fs). 

.. code-block:: bash

   # Example 1: adjust velocities for the 1st half of the velocity trajectory 
   curp conv-trj -vel \
       -p system.prmtop -pf amber \
       -i nve1.vel.nc -if netcdf --irange 1 -1 1 \
       -o stripped1.vel.nc -of netcdf --orange 1 -1 1 \
       mask -m mask.pdb

    curp conv-trj -vel \
        -p strip.prmtop -pf amber \
        -i stripped1.vel.nc -if netcdf --irange 1 -1 1 \
        -o adjusted1.vel.nc -of netcdf --orange 5 -1 5 \
        adjust-vel

Example 1 shows how to adjust the time points of the velocities to those of the atomic coordinates for the first half of the trajectory after removing unnecessary part of the system. Note that `strip.prmtop` represents the parameter/topology file generated for `mask.pdb`. As a result of the adjustment, the velocities at the time points of 10f, 20fs, :math:`\cdots`, and 4990 fs are saved to `adjusted1.vel.nc`. 

.. code-block:: bash

   # Example 2: adjust velocities for the 2nd half of the velocity trajectory 
   curp conv-trj -vel \
       -p system.prmtop -pf amber \
       -i nve1.rst -if restart --irange 1 -1 1 \
       -i nve2.vel.nc -if netcdf --irange 1 -1 1 \
       -o stripped2.vel.nc -of netcdf --orange 1 -1 1 \
       mask -m mask.pdb

   curp conv-trj -vel \
        -p strip.prmtop -pf amber \
        -i stripped2.vel.nc -if netcdf --irange 1 -1 1 \
        -o adjusted2.vel.nc -of netcdf --orange 1 -1 5 \
        adjust-vel

Similarly, example 2 shows the velocity adjustment for the 2nd half of the trajectory. Here we need to read the restart file, `nve1.rst` before reading the velocity trajectory `nve2.vel.nc`. The velocities at `t` = 4999 fs (`t` = 5001 fs) are saved in `nve1.rst` (at the 1st frame of `nve2.vel.nc`), and the velocities at `t` = 5000 fs are estimated from those at 4999 and 5001 fs. If the velocities at `t` = 5000 fs are not necessary for the final output file, you do not need to read the restart file in this example. Note that the final velocity file is generated with `"--orange 1 -1 5"` so that the first frame at `t` = 5000 fs is included in the output file named `adjusted2.vel.nc` and, thus, the velocities at 5000 fs, 5010 fs, :math:`\cdots`, and 9990 fs are save in the output file. 

Calculating interresidue heat current
======================================

Here we explain how to calculate the heat currents with

#. velocities and coordinates trajectory file
#. parameter and topology file of target system
#. configuration file for the CURP calculation

To start the calculations, please type in the following command:

.. code-block:: bash

   $ curp compute hflow.cfg > hflow.log

or 

.. code-block:: bash

   $ mpiexec -n 2 curp compute hflow.cfg > hflow.log

for parallel calculations with OpenMPI. In this case the number of cores is 2, ``hflow.cfg`` (see below)  is a configuration file for the heat current calculations and ``hflow.log`` is the log file.


These commands should produce the following two files:

*  hflow.log
*  flux_grp.nc

``flux_grp.nc`` stores the fime series of interresidue heat flux in the netcdf format.
To check the content of this file, type in the following command: 

.. code-block:: bash

   $ ncdump outdata/flux_grp.nc

Setting up ``hflow.cfg``
--------------------------

Here we show an example of ``hflow.cfg``:: 

   [input]
   format = amber
   # first_last_interval = 1 4 1
   # group_file = group.ndx

   [input_amber]
   target = trajectory
   topology_file = ../pdz3/stripped.prmtop.gz
   coordinate_format = netcdf
   coordinate_file = ../pdz3/strip.crd.nc
   velocity_format = netcdf
   velocity_file = ../pdz3/strip.vel.nc

   [curp]
   potential = amber12SB
   method = heat-flux

   group_method = residue
   flux_grain = group
   # target_atoms = 
   # enable_inverse_pair = no
   group_pair_file = gpair.ndx

   remove_trans =  yes
   remove_rotate = yes

   log_frequency = 2

   [output]
   filename = outdata/flux.nc
   format = netcdf
   decomp = no

   output_energy = no

A detailed explanation is provided below:

[input]
~~~~~~~

The input file format.

format = amber
   Read Amber formatted files.
   
first_last_interval = 1 4 1
   For the heat current calculations, the <first> and <last> frame with the interval of <intraval> frames are set in this line as <first> <last> <interval>.

group_file = group.ndx
   In this line, atom group definition file is specified. In this file, you can define an arbitrary group of atoms that is different than the standard amino acid residues.

[input_amber]
~~~~~~~~~~~~~~

In this section name, 
after ``input_`` comes the keyword specified as the format key in the ``[input]`` section.
The following keywords are used in the ``input_amber`` section.

target = trajectory | restart
   Specifies whether the input file is a trajectory file or a restart file.

topology_file = <prm_top_file>
   Set the path to the parameter and topology file.
   
coordinate_format = ascii | netcdf
   Set the format of the coordinate trajectory file.

coordinate_file = <mdcrd_file>
   File name of the coordinate trajectory file.

velocity_format = ascii | netcdf
   Set the format of the velocity trajectory file.

velocity_file = <mdvel_file>
   File name of the velocity trajectory file.

[curp]
~~~~~~~

In this section, parameters and keywords are set for the irEF calculations.

potential = amberbase | amber94 | amber96 | amber99 | amber99SB | amber03 | amber12SB
   In this line, the type of the potential function is set.

method = momentum-current | energy-flux | heat-flux
   This line specifies whether the calculation is for irEF or for atomic stress tensors (momentum current) or heat flux. In this example, we choose ``heat-flux``.

group_method = none | united | residue | file
   In this line, the unit of irEF is set.
   ``none``: No groups are defined.
   ``united``: This specifies fixed united atom groups. All of the hydrogen atoms, whether polar or apolar, belong to the united atom group represented by the heavy atom to which they attached.
   ``residue``: Groups are defined by residues unit.
   ``file``: User defined atom groups are adopted. (see ``group_file`` key in the input section.)
   none: No groups are formed.

flux_grain = atom | group | both
   Output option for the energy flow data.
   ``atom``: Output inter-atomic energy flow for all atom pairs. (not recommended) 
   ``group``: Output inter-group energy flow between all pairs of groups defined by the ``group_method`` keyword.
   ``both``: Output both of the above two data (not recommended).
   
target_atoms = 1-33
   Specifies the target segment for the calculations. In this example, atom 1 to
   33 are considered and the other atoms are neglected. If not
   specified, all atoms of the system are considered. Note that the CURP program
   excludes atoms other than the ones specifed by this option from the calculations, even when the group option is set to any of united/residue/file. 
   
group_pair_file = gpair.ndx
   Set group pair file. This option defines the set of group pair for which
   the energy flow is calculated. This can be used to focus only on
   the region of interest, saving the computational time considerably.
   Without this option, CURP calculates the energy flow between all pairs of groups.

remove_trans = yes | no
   If yes, the translational movement of the system is removed. (the default is yes)
   NOTE: Usually, we should use "yes".

remove_rotate = yes | no
   If yes, the rotational movement of the system is removed. (the default is yes)
   NOTE: Usually, we should use "yes".

log_frequency = 2
   The frequencey of output information to stdout.

[output]
~~~~~~~~

Setting the output format.

filename = outdata/flux.nc
   Filename of the heat flow data.

format = ascii | netcdf
   Format of the energy flow data. (netcdf format is highly recommended.)

decomp = no | yes
   Flag whether decompose and output the total current/flux,
   autocorrelation function and heat conductivity
   to bonded, coulomb, and van der Waals interaction.

output_energy = no | yes
   CURP is able to evaluate the energy using the atomic velocities and coordinates of the trajectory files. When set to "no", this energy value is not output. 

Calculating the time-integral of interresidue heat current
=================

After heat flux calculations, :math:`\Lambda` is calculated based on the linear response theory.
You will need the time series of :math:`\Lambda` stored in `flux_grp.nc`. Type in the following command:

.. code-block:: bash

   $ curp cal-tc \
       --frame-range 1 10 1 --average-shift 1 \
       -a outdata/acf.nc \
       -o outdata/Lambda.dat outdata/flux_grp.nc > Lambda.log

`\--frame-range <first_frame> <last_frame> <frame_interval>`
   This specifies the range of the time integration of the auto-correlation function of heat currents, :math:`\left\langle \mathbf{h}(0) \cdot \mathbf{h}(t) \right\rangle`.
   The upper and lower limits of the integral are set in <first_frame>, 
   <last_frame>, respectively. During the integration, every <frame_interval>
   frames are used for calculations.  

`\--average-shift <ave_shift>`
   In the calculation of irEC, J(0)J(t) is integrated from <first_fram> to 
   <last_frame>. Then, the origin of the integration is shifted by <ave_shift>
   and the time integration is again conducted from <first_frame> to <last_frame>.
   This procedure is then repeated until the end point of the time integration
   reaches the end of the trajectory.

`-a <file_name>`
   The time-correlation function data are output to <file_name>.
   The data format is netcdf. If this key is not specified, no data is output.

`-o <file_name>`
   :math:`\Lambda` data is output to <file_name>.

You will then obtain energy conductivity data ``output/Lambda.dat`` and the
time-correlatioin data file, ``outdata/acf.nc``.

NOTE: `\--no-axes` option should not be used here.

Format of :math:`\Lambda` data file
-------------------------

 In each line of the data file, `Lambda.dat`, a pair of residues and the corresponding value of :math:`\Lambda_{AB}` is written as <residue_A> <residue_B> :math:`\Lambda_{AB}`. The unit of :math:`\Lambda_{AB}` is measured in :math:`(\AA \rm{\cdot kcal/mol)^2/fs}`.

