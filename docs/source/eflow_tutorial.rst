=====================
Energy flow analysis
=====================

**Inter-residue energy flow (irEF) and energy conductivity(irEC)**

In 2006, we proposed the concept of inter-residue energy conductivity (irEC) to measure the strength of the residue-residue interactions (Ishikura & Yamato, *Chem. Phys. Lett.* 432: 533-537, (2006)). irEC represents the amount of the energy transferred from one residue to another per unit time and proved to be particularly useful for the detailed analysis of the intramolecular signaling mechanism of proteins. For instance, we reported the application of irEC to illustrate the photosignal transduction mechanism of a small water-soluble photosensory protein--photoactive yellow protein, as mentioned in the reference above. The value of the energy flow from :math:`j` to :math:`i` can be calculated as follows:

.. math::

   J_{i \leftarrow j} = 
    \frac{1}{2}\left( {{{\boldsymbol{v}}_i} \cdot {{\boldsymbol{F}}_{ij}}
    - {{\boldsymbol{v}}_j} \cdot {{\boldsymbol{F}}_{ji}}} \right).

The energy flow between atom groups, A and B, can be calculated by the summation of the inter-atomic energy flow between the pair:

.. math::

   J_{A \leftarrow B} =
      \sum\limits_{i \in A}^{{N_A}} {\mathop \sum \limits_{j \in B}^{{N_B}}
      {J_{i \leftarrow j}}}.

If the atom groups A and B represent residues A and B, respectively, the corresponding energy flow between A and B represents the inter-residue energy flow (irEF). Next, the energy conductivity is defined based on the Green-Kubo
linear response theory as:

.. math::

   L_{AB} = \frac{1}{{RT}}\int_0^\infty  {\left\langle
      {J_{A \leftarrow B}}(t){J_{A \leftarrow B}}(0)
      \right\rangle dt}.

In practice, the time integral is replaced by the sum of the auto-correlation function of :math:`J` at discrete time points. We came to know, from our experience, that the time resolution of :math:`J` should be no larger than 10 fs for the accurate estimation of the irEC. 
The CURP program performs these calculations for the proteins
and visualizes the energy conductivities using an energy exchange network (EEN).This allows the user to grasp the pattern of the inter-residue interactions at a glance. In principle, it is possible to calculate irECs based on a single NVE trajectory. However, it is highly recommended to calulate thermal averages of irECs based on multiple NVE trajectories that start from different points on the energy landscape of the target protein. If the number of the NVE trajectories are sufficiently large, the thermally averaged irECs are fair representations of the irECs of the thermally fluctuating protein system.  

.. note:: In the previous paper(*CPL* ,2006), we calculated the irEC based on the atomic decomposition of the potential energy using classical force-field functions, where the way of the decomposition involves arbitrariness. In the CURP program, the irEC calculation is performed in a different manner based on the pairwise decomposition of forces (*CPL* 539(2012)144). Accordingly, the numerical results of the irEC may not be coincide with those obtained earlier.

We will illustrate how to perform energy flow analysis using the PDZ3 domain as an example. To calculate the irEC described above, coordinates and velocities trajectory of the system of interest are required. The sample trajectory data of the PDZ3 domain can be downloaded from the following link:

`PDZ3 domain trajectories <http://www.comp-biophys.com/resources/md-trj/pdz3.tar.gz>`_

Note that each trajectory file can be larege (a few GB).
Below is a summary of the MD simulation:

*  We used the Amber12 program.
*  Multiple snapshots were extracted from the NPT trajectory and atomic cooridinates and velocities were saved for restart conditions, each of which was used to conduct the subsequent NVE simulation for 1 ns.
*  During the NVE simulation, atomic velocities and coordinates were saved every 10fs, thus the total # of frames was (1ns/10fs =) 100,000.
*  Since we are interested in the inter-residue energy flow, we removed the solvent degrees of freedom from the trajectory files.

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

`--irange`
   Specifies <first_frame>, <last_frame> and <frame_interval> for the input trajectory file.
   `-1` for <last_frame> represents the last frame of the trajectory.
   For example, if "1 -1 5" is provided, the conv-trj command reads
   1st, 6th, 11th, ... , up to the last frames out from the trajectory.

`-o`
   Specifies the trajectory file for output. You are not allowed to
   specify this option multiple times.

`-of`
   Specifies the format of the output trajectory file.

`--orange`
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

Calculating irEF
=================

Here we explain how to calculate irEF using PDZ3 as an example, with

#. velocities and coordinates trajectory file (see above)

   *  You will find test data in `tutorial/pdz3-eflow/amber-mddata`.
   *  Atomic coordinates and velocities saved every 10 fs and the total number frames is 100.

#. parameter and topology file of target system
#. configuration file for the CURP calculation

To start the calculations, please type in the following command:

.. code-block:: bash

   $ curp compute eflow.cfg > eflow.log

or 

.. code-block:: bash

   $ mpiexec -n 2 curp compute eflow.cfg > eflow.log

for parallel calculations with OpenMPI. In this case the number of cores is 2, ``eflow.cfg`` (see below)  is a configuration file for the irEF calculations and ``eflow.log`` is the log file.

Alternatively, ``run_eflow.sh`` performs the equivalent tasks.

.. code-block:: bash

   $ cd tutorial/pdz3-eflow/eflow+ec
   $ ./run_eflow.sh

These commands should produce the following two files:

*  eflow.log
*  flux_grp.nc

``flux_grp.nc`` stores the fime series of irEF in the netcdf format.
To check the content of this file, type in the following command: 

.. code-block:: bash

   $ ncdump outdata/flux_grp.nc

Setting up ``eflow.cfg``
--------------------------

Here we show an example of ``eflow.cfg``:: 

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
   method = energy-flux

   group_method = residue
   flux_grain = group
   # target_atoms = 
   # enable_inverse_pair = no
   group_pair_file = gpair.ndx

   remove_trans =  no
   remove_rotate = no

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
   For the irEF calculations, the <first> and <last> frame with the interval of <intraval> frames are set in this line as <first> <last> <interval>.

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

method = momentum-current | energy-flux
   This line specifies whether the calculation is for irEF or for atomic stress tensors (momentum current). In this example, we choose ``energy-flux``.

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
   If yes, the translational movement of the system is removed.

remove_rotate = yes | no
   If yes, the rotational movement of the system is removed.

log_frequency = 2
   The frequencey of output information to stdout.

[output]
~~~~~~~~

Setting the output format.

filename = outdata/flux.nc
   Filename of the energy flow data.

format = ascii | netcdf
   Format of the energy flow data. (netcdf format is highly recommended.)

decomp = no | yes
   During the calculations, choose whether the energy is decomposed into different
   components.

output_energy = no | yes
   CURP is able to evaluate the energy using the atomic velocities and coordinates of the trajectory files. When set to "no", this energy value is not output. 

The log file looks like  `this <./curplog.txt>`__.

Calculating irEC
=================

After irEF calculations, irEC is calculated based on the linear response theory.
You will need the time series of irEF stored in `flux_grp.nc`. Type in the following command:

.. code-block:: bash

   $ curp cal-tc \
       --frame-range 1 10 1 --average-shift 1 \
       -a outdata/acf.nc \
       -o outdata/ec.dat outdata/flux_grp.nc > ec.log

`--frame-range <first_frame> <last_frame> <frame_interval>`
   This specifies the range of the time integration of the auto-correlation function of irEF, :math:`(J(0)J(t))`.
   The upper and lower limits of the integral are set in <first_frame>, 
   <last_frame>, respectively. During the integration, every <frame_interval>
   frames are used for calculations.  

`--average-shift <ave_shift>`
   In the calculation of irEC, J(0)J(t) is integrated from <first_fram> to 
   <last_frame>. Then, the origin of the integration is shifted by <ave_shift>
   and the time integration is again conducted from <first_frame> to <last_frame>.
   This procedure is then repeated until the end point of the time integration
   reaches the end of the trajectory.

`-a <file_name>`
   The time-correlation function data are output to <file_name>.
   The data format is netcdf. If this key is not specified, no data is output.

`-o <file_name>`
   Energy conductivity data is output to <file_name>.

A useful script file, ``run_ec.sh`` is available for these calculations:

.. code-block:: bash

   $ cd tutorial/pdz3-eflow/eflow+ec
   $ ./run_ec.sh

You will then obtain energy conductivity data ``output/ec.dat`` and the
time-correlatioin data file, ``outdata/acf.nc``.

Format of irEC data file
-------------------------

 In each line of the data file, `ec.dat`, a pair of residues and the corresponding value of irEC is written as <residue_A> <residue_B> <L_AB * RT>, where <A> = donor residue, <B> = acceptor residue, L_AB = irEC between the residues A and B, R is the gas constant, and T is the absolute temperature (= 300 K for most cases). The unit of <L_AB * RT> is measured in :math:`(kcal/mol)^2/fs`. Note that the order of <A> and <B> makes no difference bacause the value of L_AB is evaluated for the pair (A,B) without any directionality. However, the difference between the donor and acceptor could be important in some cases. For example, the interatomic electron tunneling current has directionality and in that case the order of the donor and the acceptor is important.

Drawing the EEN(Energy Exchange Network)
=========================================

Example scripts to draw the EEN are found in ``tutorial/pdz3-eflow/een`` directory, that contains:

graph-een
   This directory includes some useful scripts that: (1) remove the residue pairs adjacent in the sequence, (2) renumber the residue numbers.

graph-een-apo
   Draw an EEN graph for apo PDZ3 domain.

graph-een-a3rem
   Draw an EEN graph for :math:`\alpha 3`-truncated PDZ3 domain.

graph-een-diff
   Draw a difference graph of the EEN graphs of apo- and :math:`\alpha 3`-truncated PDZ3 domains.

view-een-3D
   Mapping EEN connectivity on the 3D structure using PyMOL.

Drawing the EEN
----------------


Preparation
~~~~~~~~~~~

In this example, the working directory for drawing EEN graph is `graph-een`. We used the data file `apo.ec.dat` in this directory. If user wish to use the user's own irEC data file, perform the following steps: 

#. copy the sample graph-een directory to a location wherever you like.
#. copy the user's irEC data file to the new graph-een directory.
#. Set the `ec_fp` variable in the `config.mk` file to the user's irEC data file name.

Then you will get the EEN graph output by running `make`.


Basic Usage
~~~~~~~~~~~

After editing the `config.mk` file and specifying some parameters for the
`graph-een` command, run `make` to obtain an `strong.ec.pdf` file, which 
graphically illustrates the EEN.

Other Usage
~~~~~~~~~~~~

By running make with different options, you can obtain the
EEN in different representations as follows:

`make`
   The standard way to generate the EEN of the target system.
   Both weak interactions and strong interactions
   are shown on the graph.

`make clean`
   Clean up the working directory. EEN graph output files
   and the associated intermediate files are removed.

`make strong`
   Only strong inter-residue interactions that are greater than
   the thresh value (see below) are shown on the EEN graph.
   
`make weak`
   Weak interactions that are greater than the thresh_weak value (see below)
   are shown on the graph together with the strong interactions.

Setting config.mk
~~~~~~~~~~~~~~~~~~~

To customize the graph drawing conditions, the `config.mk` fle should be edited.
An example of the `config.mk` file
is shown below: 

ec_fp = ../all.ec.dat
   The name of the file that contains the irEC data. 

fix_resnums = 1:306
   Renumbering of the residue numbers of the polypeptide chain. In this example, 
   the number of residue 1 is changed to 306. You may add multiple items
   separated by a space between them.

thresh = 0.008
   Setting the threshold value for drawing the EEN graph with the
   strong option (see above).  

thresh_weak = 0.003
   Setting the threshold value for drawing the EEN graph with the
   weak option (see above).  

line_values = 0.015  0.008  0.003
   The threshold values for line attributes. Number of elements in the list must be equal with all line attribute.

line_colors = red  blue  green
   The colors of lines. Number of elements in the list must be equal with all line attribute.

line_thicks = 4.0  4.0  2.5
   The thickness of lines. Number of elements in the list must be equal with all line attribute.

line_weights = 5.0  3.0  1.0
   The weight of line. Number of elements in the list must be equal with all line attribute.
   
other_opts = --ratio 0.3 --direction TB -I
   Setting other options passed to the `graph-een` command. (see below)

Selecting the important residue pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. HERE

  The CURP program calculates irEC of a target protein. The `all.ec.dat` files stores the data, which is then processed by other scripts, such as `renum-residue`, `sele noasa`, ... from `curp tool` command, selecting the irEC on which the user wants to focus.

.. code-block:: bash

   $ cat ../all.ec.dat | ./renum_residue.py 1:10

In this directory, a useful make target, `make sel.ec.dat`, is provided to run multiple python scripts at a time.

.. .. code-block:: make

.. **Data files**
.. 
.. - asa.dat
..    This file stores the list of exposed residues with solvent accessibility greater than 0.3.
.. 
.. - ss-b2AR.dat
..    This file describes the secondary structures of the target protein.
.. 
.. - label_conv.dat
..    This file describes the renaming rules of amino acid residuesor hetero groups. One rule is given per one line. ``#`` indicates comment.

Brief usages of scripts that process irEC data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``sed -e "s/CYX/CYS/g" < ec-org.dat > ec-new.dat``
   Convert the cysteine name from CYX used in Amber into the standard one.

.. sel_thresh_dist.py dist.dat 6 3 < ec-org.dat > ec-new.dat
   Select only irECs the distance of which have less than 6.0 Å.
   The second argument means that the third column is applied
   as thereshold value.

``curp tool renum-residue $(fix_resnums) < ec-org.dat > ec-new.dat``
   Renumbering the residue numbers according to the `fix_resnums` variable.

``curp tool noneighbor WAT $(ligand) < ec-org.dat > ec-new.dat``
   Remove neighboring residue pairs that indicate covalent peptide bonds.

``curp tool nocap WAT $(ligand) < ec-org.dat > ec-new.dat``
   Remove residue pairs that contain the capped residues.

.. - convert_labels.py
..    第一引数で指定したファイルの中で定義された残基名を変換することが出来る。
..    元は `graph_ec.py` に含まれていたが、
..    他のスクリプトでも使えるようにするために、
..    単純なスクリプトとして分離した。
.. 
.. - select_noasa.py  asa.dat  <  ec-org.dat  >  ec-new.dat
..    Remove the exposed residues listed in `asa.dat` from the original data of irEC, `ec-org.dat`.
.. 
.. - select_nohelix.py  ss-b2AR.dat  13  <  ec-org.dat > ec-new.dat
..    Remove the hydogen bonded residue pairs of (i,i+4)-type in helices. The first argument specifies the file name that describesthe secondary structure, while the second argument describes the shift of residue number. In this case, residue i in `ec-org.dat` and `ec-new.dat` corresponds to residue i+13 in ss-b2AR.dat.
.. 
.. - select_noneighbor.py  WAT CAU < ec-org.dat > ec-new.dat
..    Remove the neighboring pairs of residues. The exceptions of this rule apply to the residues whose names are specified in the second, third, ... arguments.
.. 
.. - convert_names.py label_conv.dat < ec-org.dat > ec-new.dat
..    Residue name conversion is performed. The rules are described in `label_conv.dat`.
.. 
.. - renum_residue.py 1:14 100:150  ... < ec-org.dat > ec-new.dat
..    残基番号を変更する。
..    この例では次のように変更される:
.. 
..       | 1 => 14
..       | 2 => 15
..       | .
..       | .
..       | 99 => 112
..       | 100 => 150
..       | 101 => 151
..       | .
..       | .
.. 
.. convert_BW.py  b2AR-BW.dat < ec-org.dat > ec-new.dat
..    B&W表記のデータベースを用いて、データをB&W表記に変更する。

Usage of `graph-een`
~~~~~~~~~~~~~~~~~~~~~~~

.. HERE

Here we describe parameters for the ``curp graph-een`` command.

`-h, --help`
   Display help menu.

`-f, --output-een-filename`
   Specifies the output file name for the EEN.
   An obtained figure in pdf format can be edited with Adobe Illustrator.

`-r, --threshold`
   Specifies the threshold value of irEC for drawing. If a residue pair (A, B) has energy conductivity L_AB such that L_AB * RT is greater than this threshold value, then the residue pair appears in the EEN. Here R is the gas constant and T is the absolute temperature (= 300K for most cases). The default value is 0.01 :math:`(kcal/mol)^2/fs`.

`--ratio`
   Specifies the inverse aspect ratio. This value specifies the ratio of the height of the image to its width. Without this parameter, the inverse aspect ratio is set automatically.

`-c, --cluster-filename`
   The cluster structure of the EEN is illustrated in the image. The cluster structure of the nodes (residues) in the EEN is described in the file specified by --cluster-filename. See `Format of cluster definition file` for details.

`-n, --node-style-filename`
   Specifies the file name for the node style definition. See `Format of the node style definition file` for details.

`--direction`
   Specifies the orientation of the EEN, which is drawn from left to right with **LR** and top to bottom with **TB**. The default value is **TB**.

.. -t, --targets
..    Specify the target residues for drawing the image of the EEN.
..    This parameter describes the residue numbers for drawing the image of the EEN in either of the following ways.
..    
..    *  1-
..    *  -340
..    *  5-8
..    *  -9 22 50 80 93-108

The selected residues will appear in the image irrespective of whether the residues are donors or acceptors.

`--with-one-lettr`
   Specfication of this option leads to the showing of one letter symbols  
   for the amino acid residues.

.. --forced-output-nodes
..    The communication map forcelly includes given residue number's nodes even if
..    their nodes have the irEC values that are less than threshold value.


`-p,  --bring-node-pair-close-together`
   Bring the node pair close together on the EEN graph. ex), 1:2, 3:5, 3:15.

.. `-I, --hide-isolated-nodes`
   .. Hide isolated nodes when applying multiple EC files.



Input Data
~~~~~~~~~~~

node.cfg
   This file describes the style of each node in the image.
   See `Format of node style definition file`_ for details.

cluster.cfg
   This file describes the cluster structure of the nodes in the EEN. See `Format of the cluster definition file`_ for details.
    
Format of the cluster definition file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we describe the format of the cluster definition file, `cluster.cfg`.

Example of `cluster.cfg`::

   [Ligand binding pocket]
   322 323 324 325 326 327 328 331 339 372 376 379 380

   shape = box, filled
   color = black
   fillcolor = yellow
   fontsize = 30
   .
   .
   .

In the above example each cluster name is indicated by square brackets ``[`` and ``]``. This is followed by the listing of the members of each cluster. For instance, residues 29-60 are included in the cluster TM1.

Cluster attributes may be specified as well. You can specify multiple attributes 
such as `color` (= red, yellow, ...), `shape` (= box, rounded,...) and so on. Basically, any `Graphviz-attributes`_ can be specified.

Examples of common attributes:

*  style =
*  shape = box, rounded, circle, filled, ...
*  Make a backgroud white

   *  color = white
   *  fillcolor = white

*  fontcolor = -


.. note::
   Sometimes, you may group multiple nodes together. To do so::

      fillcolor = white
      label =  
      color = white

   An important point is to leave `label`-attribute blank.

.. _Graphviz-attributes: http://www.graphviz.org/doc/info/attrs.html


Format of node style definition file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file describes node attributes. The format of this file is similar to that of cluster definition file. The difference between the two files is that the cluster definition file describes the attributes of the cluster as a whole, while the node style definition file specifies the attributes of each individual node.

Example of the node style definition file::

   [default]
   1-

   fillcolor = lavender
   fontcolor = black
   fontname = Arial Bold
   # fontname = Helvetica Neue Bold
   fontsize = 14
   height = 0.3
   width  = 0.7

   [alpha helix 3]
   393-399

   fillcolor = black
   fontcolor = white

A section name is indicated by square brackets ``[`` and ``]``. In contrast to the cluster names, section names do not actually appear on the image of the EEN. If a node belongs to two different sections, the attribute setting of the new section overwrites the previous one. Therefore, it is convenient to put the ``[default]`` section at the top and set the default attributes there. The line, ``1-``, in the default section means that the default attributes apply to all of the nodes.

Attribute examples:

*  label = <value>
*  style = <value>
*  shape = elliple | box | rounded | circle |
*  color = <value>
*  fillcolor = <value>
*  fontcolor = <value>
*  fontsize = <value>
*  height = <value>
*  width = <value>

Drawing the differential EEN
-----------------------------

The differential EEN graph illustrates the difference of the irEC between a pair of homologous proteins or different states of the same protein as an EEN representation. The working directory for this is `graph-een-diff`. In this example, we show how to draw the differential EEN between the apo- and the :math:`\alpha 3`-truncated forms of PDZ3 domain. The respective irEC data are stored in `chart-een-apo/outdata/sel.ec.dat` and `chart-een-a3rem/outdata/sel.ec.dat`. Before drawing the differential EEN graph, set the `ec_fp1` (`ec_fp2`) variable to `chart-een-apo/outata/sel.ec.dat` (`chart-een-a3rem/outdata/sel.ec.dat`) in the `config.mk` file.

The differential EEN is consisted of two parts: For the 1st (2nd) part, the values of irEC increase (decrease) from `ec_fp1` to `ec_fp2`, and these parts are denoted as ``inc`` and ``dec``, respectively. Note that ``inc`` (``dec``) would contain the node pairs that are found only one of either `ec_fp1` or `ec_fp2`. In such cases, there are two kinds of possibilities as described below:

Union 
   We consider such pairs for the differential EEN. The missing pairs in either `ec_fp1` or `ec_fp2` are assumed to have the virtual value (see below) specified by the `dval` variable in `Makefile`. Such pairs are shown in the `inc-fill.ec.pdf` (`dec-fill.ec.pdf`) file.

Intersection
   Such pairs are neglected for the differential EEN.

In `Makefile` in the `graph-een-diff` directory, you will the following variables:

dval = 0.0001

Set the virtual value of irEC if a node pair is missing in either `ec_fp1` or `ec_fp2` (see above). If you running `make`, you will get the following five files:

`inc.ec.pdf`
   The `inc` part of the differential EEN (Union).

`dec.ec.pdf`
   The `dec` part of the differential EEN (Union).

`inc-fill.ec.pdf`
   The `inc` part of the differential EEN (Intersection).

`dec-fill.ec.pdf`
   The `dec` part of the differential EEN (Intersection).

`both.ec.pdf`
   Both `inc` and `dec` parts are shown as solid and dashed lines, respectively. (Union)

Sometimes, you may be interested in only on of the five files.
In such a case, for instance, if you type the following command:

.. code-block:: bash

   $ make inc

then you will get only `inc.ec.pdf`.


.. EEN and 3D structure
.. ----------------------
 
.. .. todo:: This section should be written.
.. 
.. **Directory** : view-een-3D
.. 
.. It is possible to chart the EEN on the 3D structures with
.. the PyMOL graphics program. Here we provided some useful python scripts in
.. this directory.
.. Note that these scripts need the PyMOL program installed in your local environment.
.. 
.. pml file
.. --------
.. You need an appropriate pml file to specify the setup environment for the PyMOL program.
.. 
.. Usage
.. -----
.. 
.. It is possible to chart the EEN on the 3D structure by running the following command.
.. 
.. .. code-block:: bash
.. 
..    $ pymol ./system.pdb b2ar_ec.pml -r view_ec.py -- ec-sel-proper.dat rank.dat
.. 
.. In this example, `system.pdb` and `b2ar_ec.pml` represents the PDB file of the protein and the setup file of the PyMol session. The `view_ec.py` is the main part of the python script, and this script reads two files `ec-sel-proper.dat` and `rank.dat`, where the former stores the energy conductivity data and the latter stores the evolutionary trace data, respectively. For convenience, this command line can be executed by simply running the `make` command after editing `config.mk` that controls appearance behavior.
.. 
.. .. code-block:: bash
.. 
..    $ make
..    
..       or
.. 
..    $ make view
.. 
.. Moreover, it is possible to create the movie based on the 
.. 
.. .. code-block:: bash
.. 
..    $ make movie
.. 
.. To clean up
.. 
.. .. code-block:: bash
.. 
..    $ make clean
.. 
.. .. note::
..    The residue numbering scheme of the donor and acceptor in the energy
..    conductivity data file, `ec-sel-proper.dat` should be coinside with that
..    in the PDB file, `system.pdb`. The `data-ec/convert-names.py` script may be
..    helpful for editing the residue numbers for this purpose.
.. 
.. .. ここで使われるエネルギー伝導度のデータは、
..    最初に読み込まれる構造ファイルの残基番号に依存している。
..    そのため、構造ファイル内での残基番号に合致するように
..    伝導度データのファイルを作成すること。
..    これには、 `data-ec/convert-names.py` を上手く使って、
..    `view_ec.py` のためのエネルギー伝導度データを特別に作ると良い。
.. 
.. .. todo:: スクリプトを変更して、残基番号に気を付ける文書を改正する。

