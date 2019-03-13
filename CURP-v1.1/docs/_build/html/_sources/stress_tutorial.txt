=============================
Atomic stress tensor analysis 
=============================

CURP performs atomistic stress tensor analysis. According to Hardy's forluma,
the Cauchy stress tensor of the *i*\ th atom is given by:

.. .. math::

   \boldsymbol{\sigma}_i = \frac{1}{\Omega_i} \left(
   \frac{1}{2} \sum_j \boldsymbol{F}_{ij} \otimes \boldsymbol{r}_{ij}
   + m_i \boldsymbol{v}_i \otimes \boldsymbol{r}_{ij}
   \right)

.. image:: stress.*
   :align: center
   :height: 60px

where :math:`\Omega_i` is the volume of atom :math:`i`,
:math:`\boldsymbol{F}_{ij}` is the pairwise interatomic force
acting on :math:`i` from :math:`j`,
:math:`\boldsymbol{r}_{ij}` is the relative displacement vector
from atom :math:`j` to atom :math:`i`,
:math:`m_i`  is the mass of atom :math:`i`, and
:math:`\boldsymbol{v}_i` is the velocity vector of atom :math:`i`.

As can be seen from the above formula, the atomic stress tensor of atom :math:`i` depends on both the
position and the velocity trajectories during the MD simulations.

The calculation procedure consists of the following steps:

.. contents::
   :local:
   :depth: 1

Obtaining the MD trajectory
===========================

We show two examples of calculations for **ALA-ALA-ALA in vacuum/solvent**.
NVT MD simulations were performed using the Amber12 program with the time step
of 0.5 fs at room temperature. In addition to the trajectory files,
a parameter topology file is required for each system.

Stress tensor calculation
=========================

Examples for the alanine tripeptide in the vacuum and solvent environments
are found in the directories ``ala3invac-stress`` and ``ala3insol-stress``,
respectively.
Each directory contains three subdirectories named (1) united, (2) residue and (3) whole.

In practice, it is possible to calculate atomic stress tensors for all atoms including the hydrogen atoms, and CURP considiers all atoms in the calculations. However, the all-atom representation is often inconvenient for the stress tensor analysis because the values of the atomic volume obtained by the Voronoi tessellation are sometimes quite small for the hydrogen atoms, leading to large and significant fluctuations.
To overcome this issue, a united-atom approach has been employed. 

.. math::
   
   \boldsymbol{\sigma}_{J} = 
      \frac{1}{{{\Omega _J}}}\sum\limits_{i \in J}
      {\left( {\frac{1}{2}\sum\limits_{k \ne i}
      {{{\boldsymbol{F}}_{{\boldsymbol{ik}}}} \otimes {{\boldsymbol{r}}_{{\boldsymbol{ik}}}}} }
      \right)}

,where :math:`J` represents a group of atoms. Here, each unit of the atom group could be either (1) united atom (heavy atom together with the attached hydrogen) or (2) amino acid residue, or (3) the whole molecule.

Running the program
-------------------

First, change the directory to the working directory for the stress tensor analysis.

.. code-block:: bash

   $ cd $CURP_HOME/tutorial/ala3invac-stress/united

To start the calculations, type the following command:

.. code-block:: bash

   $ $CURP_HOME/bin/curp < run.cfg > log

where ``run.cfg`` (see below for details) is the configuration file 
for the analysis. An alternative way is to run the following script:

.. code-block:: bash

   $ ./run.sh

The calculations will be completed shortly and the prompt will appear at the shell.
The results are stored in the files in two forms: (1) ``stress_atm_dat00000``:
individual stress tensors for every atom, (2) ``stress_grp_dat00000``:
averaged stress tensors for each group, where a group can be either united atom, residue, or the whole molecule.

.. note::
   The following lines appear at the bottom of the log file.
   **Detailed timing in calculator object**,
   **Detailed timing for parallel prossing**,
   **The summary of the elapsed times** 
   The timging information is currently incorrect. Please ignore these lines for now.

Setting run.cfg
---------------

An example of `run.cfg` is shown below.

::

   [input]
   format = amber
   first_last_interval = 1 200 1
   # group_file = group.ndx

   [input_amber]
   target = trajectory
   topology_file = ../system.prmtop
   coordinate_format = ascii
   coordinate_file = ../sam.mdcrd.gz
   velocity_format = ascii
   velocity_file = ../sam.mdvel.gz

   [volume]
   method = voronoi
   output_volume_file  = outdata/volumes.dat
   output_gvolume_file = outdata/gvolumes.dat
   voronoi_cutoff = 6.0
   voronoi_no_hydrogen = yes
   voronoi_solvation = RANDOM20
   voronoi_probe_length = 2.8
   # voronoi_output_solvation_file = outdata/solvated.pdb

   [curp]
   potential = amber99
   method = momentum-current

   group_method = united
   decomp_group_current = no
   target_atoms = 1-33

   remove_trans = yes
   remove_rotate = yes
   log_frequency = 10

   [output]
   filename = outdata/stress.dat
   decomp = yes
   frequency = 10000


[input]
~~~~~~~

In this section, we specify the input file format.

format = amber
   CURP reads the input files in the Amber format.
   
first_last_interval = 1 200 1
   Three numbers, <first>, <last>, and <interval> are given in this order, to specify that CURP will read the MD frames from the <first> to the <last> position at intervals of <interval> frames. 

group_file = group.ndx
   Specifies the name of the 'group definition file'. Each group is defined as a list of atoms. A group could be either a residue or an arbitrary set of atoms.

[input_amber]
~~~~~~~~~~~~~

This section name depends on the `format` keyword specified in the ``[input]`` section. The name of this section is `[input_amber]` because the `format` keyword in the ``[input]`` section is "amber".

For the ``[input_amber]`` section, the following parameters must be specified.

target = trajectory | restart
   If `trajectory` is specified, CURP assumes that the input files are trajectory files.

topology_file = ../system.prmtop
   Parameter and topology file.

coordinate_format = ascii | netcdf
   Format of the trajectory files.

coordinate_file = <mdcrd_file>
   The name of the coordinates trajectory file.

velocity_format = ascii | netcdf
   Format of the velocity trajectory files.

velocity_file = <mdvel_file>
   The name of the velocity trajectory file.

[curp]
~~~~~~

In this section, we set the parameters for the calculations of pairwise interatomic forces.

potential = amberbase | amber94 | amber96 | amber99 | amber99SB | amber03 | amber12SB
   Type of force-field functions.

method = momentum-current | energy-flux
   `momentum-current`: Atomic stress tensors analysis.
   `energy-flux`: Interatomic energy flow analysis.

group_method = none | united | residue | file
   Unit of atom groups.
   ``none``: No groups are defined.
   ``united``: Each heavy atom, either polar or nonpolar, represents a united atom "group" and all hydrogen atoms attached to the heavy atom  belong to the "group".
   ``residue``: Groups are defined by residues unit.
   ``file``: Groups are defined by the file specified by the ``group_file`` keyword in the input section.

target_atoms = 1-33
   The target region for the analysis.
   In this case, atomic stress tensors are calculated for the 1st to the 33rd atoms.

remove_trans = yes | no
   Translational motions are removed from the trajectory, if this keyword is set to "yes" (strongly recommended).

remove_rotate = yes | no
   Rotational motions are removed form the trajectory, if this keyword is set to "yes" (strongly recommended).

log_frequency = 10
   Print results to stdout every log_frequency steps.

[volume]
~~~~~~~~

For the atomic stress tensor analysis, CURP performs Voronoi tessellation
to calculate the atomic volumes.

method = voronoi | none | vdw | outer
   ``voronoi``: Atomic volume calculations are performed by Voronoi tesselation.
   See `Specification of curp setting file` in the reference section for other values.
   We coded a Fortran90 program for the atomic volume calculation based on Voronoi 
   tesselation algorithm found in the literature ("Computer Simulation of Liquids", 
   M.P. Allen and D.J. Tildesley eds., Oxford Univ. Press (1987)).

output_volume_file = <file_path>
   The path to the output file for the atomic volumes.

output_gvolume_file = <file_path>
   The path to the output file for the group volumes.

voronoi_cutoff = 6.0
   The cutoff length for the neighbor candidate search in voronoi tessellation.

voronoi_no_hydrogen = yes | no
   Treatment of the hydrogen atoms in Voronoi tessellation.
   If this option is "yes", then the hydrogen atoms are neglected for the Voronoi tessellation and the evaluation of the group stress tensors, even when a hydrogen atom is included in the target for the calculation in the [curp] section. Within the program, in this case, the atomic volume of a hydrogen atom is set to be 8.0 Å :sup:`3` to avoid division by zero for convenience.

voronoi_solvation = RANDOM20 | none
   For atoms exposed on the surface of a target protein molecule
   in vacuum, it is not possible to define the Volonoi polyhedron. 

   `RANDOM20`: To avoid this problem, the protein molecule is placed
   in a hydration sphere (radius = 20 Å) of randomly generated 
   non-overlapping water molecules. If your target system is so large
   that the hydration sphere is not able to cover the system, plase use
   the `none` option and solvate the system with a sufficiently large number of 
   hydration layers. (see below)

   `none`: This keyword is specified when CURP performes the stress
   tensor analysis for the system including explicit solvent layers.

voronoi_probe_length = 2.8
   The probe length [Å] of the solvation for the Voronoi method.
   If voronoi_solvation is set, CURP removes water molecules within
   <voronoi_prove_length> [Å] from all atom of the system.

[output]
~~~~~~~~

Output format specification.

filename = outdata/stress.dat
   Specify the naming convention of the datafile for the stress tensor analysis.
   In this example, the file name is `outdata/stress_grp.dat00001`,
   `outdata/stress_grp.dat00002` , ... and so on.

   Please note that the datafile is generated for each force component if `decomp = yes`. (see below)

frequency = 10000
   The number of frames output to a single datafile. It is highly recommended to estimate the file size of a single datafile in advance.

decomp = no | yes
   If "yes", stress tensors are separated into different components.

Setting group.ndx
-----------------

As mentioned above, the user can define any group of atoms for the
stress tensor analysis. When the `group_method` keyword is set to be
`file` in the [curp] section, and `group_file` name is specified in the
[input] section, user defined grouping is applied.

In the following example, stress tensor analysis is performed separately
for the main chain parts and for the side chain parts.

::

   [01_ALA_M]
   1-6 11-12

   [01_ALA_S]
   7-10

   [02_ALA_M]
   13-16 21-22

   [02_ALA_S]
   17-20

   [03_ALA_M]
   23-26 31-33

   [03_ALA_S]
   27-30

Here, each group name is surrounded by `[` and `]`.
Then the list of the individual atoms is given in which each atom is separated
by a space, tab or a blank line. Alternatively, a range of atoms can
be indicated by using `-`.

Output files
------------

stress tensor data file
~~~~~~~~~~~~~~~~~~~~~~~

For each snapshot of the MD trajectory, 3x3 stress tensors are output to
`outdata/stress_grp.dat00000`. An example is shown below:

.. image:: stress_grp.png
   :align: center

.. .. code-block:: bash

..   $ cat outdata/stress_grp.dat00000
   %title  
   %format time {:>15.3f} [ps]
   %format data {name:>14s}  {xx:012.7e} {xy:012.7e} {xz:012.7e} {yx:012.7e} {yy:01
   2.7e} {yz:012.7e} {zx:012.7e} {zy:012.7e} {zz:012.7e} 
   %time          0.000 [ps]
   %data
          00001_N  1.9281438e-01 -1.2271606e+00 6.1073072e-01 -1.2271606e+00 2.9022 752e-01 -6.5915250e-01 6.1073072e-01 -6.5915250e-01 5.8197837e-01 
         00005_CA  8.6272863e-01 -3.5911721e+00 -5.7126894e-01 -3.5911721e+00 1.586 1648e-01 4.6082898e-01 -5.7126894e-01 4.6082898e-01 -2.8353371e-01 
         00007_CB  -1.9837408e+00 6.3024438e-01 1.2160073e-01 6.3024438e-01 -5.5510 539e-01 -4.1604746e-01 1.2160073e-01 -4.1604746e-01 6.1241842e-02 
          00011_C  -1.0033768e+01 5.0524848e+00 -7.7353054e-02 5.0524848e+00 -3.054 9891e+00 2.3283834e+00 -7.7353054e-02 2.3283834e+00 -1.5796993e+01 
          00012_O  1.0840796e+00 -5.1829702e-01 -3.0867468e-01 -5.1829702e-01 3.883 4391e-01 -4.1876701e-01 -3.0867468e-01 -4.1876701e-01 3.0278155e+00 
   .
   .
   .
   %time          1.000 [ps]
   %data
   .
   .
   .

Comment lines begin with ``%``.
For each time frame, nine elements (xx, xy, xz, yx, yy, yz, zx, zy, zz) are printed for each group.

voluemes.dat, gvoluems.dat
~~~~~~~~~~~~~~~~~~~~~~~~~~

Atomic volumes and group volumes are printed to the files
specified by the `output_volume_file` and
`output_gvolume_file` keywards in the [volume] section.
If no file names are specified, no data is printed.

Analysis of stress tensor
=========================

The direct output of the stress tensor analysis itself is inconvenient for 
two reasons: (1) the output file contains a large amount of data 
of the time-series of the MD trajectory and (2) the stress tensor has 3x3=9
components.

To simplify the analysis,
we show how to diagonalize and time average the stress tensor as follows:

Time-averaging:

.. math::
   \boldsymbol{\sigma}_{J} = \left\langle
      {\frac{1}{{{\Omega _J}}}\sum\limits_{i \in J} {\left( {\frac{1}{2}\sum\limits_{k \ne i} {{{\boldsymbol{F}}_{{\boldsymbol{ik}}}} \otimes {{\boldsymbol{r}}_{{\boldsymbol{ik}}}}} }
      \right)} } \right\rangle

Diagonalization:

.. math::
   \boldsymbol{\sigma '} = \begin{pmatrix}
      \sigma _{x'x'}&0&0 \\
      0&\sigma _{y'y'}&0 \\
       0&0&\sigma _{z'z'}
   \end{pmatrix}

Root-mean-square amplitude of the diagonal components:

.. math::
   \bar{\sigma '} =
      \sqrt {\sigma _{x'x'}^2 + \sigma _{y'y'}^2 + \sigma _{z'z'}^2}

The last quantity is used to measure the magnitude of the stress tensor.

These calculations are performed by the `ana.sh` script.

.. code-block:: bash

   $ ./ana.sh

`ana.sh` uses  `$CURP_HOME/script/simplify_tensor.py`,
which is executed as an argument to `$CURP_HOME/bin/ana-curp`:

.. code-block:: bash

   $ $CURP_HOME/bin/ana-curp simplify_tensor.py


Result
-------

Finally, a detailed summary of the calculated data is output to `grp-sim.dat`.

.. image:: ./grp-sim.png
   :align: center

..   #label id  name        total      bond      angle    torsion   improper    coulomb        vdw    kinetic
       00001  N       0.6027320      0.29       0.61       0.05       0.00       0.95       1.07       0.10
       00005  CA      4.4768001      0.64       1.42       1.79       3.37       0.39       0.04       0.12
       00007  CB      0.1212715      0.23       0.14       0.05       0.00       0.07       0.03       0.09
       00011  C      20.1206626      0.56       1.11       3.54      17.40       2.94       0.05       0.10
       00012  O       3.2316832      0.11       0.50       1.29       2.18       0.33       0.07       0.03
       00013  N       4.2103065      1.01       0.96       2.38       3.36       0.75       1.00       0.12
       00015  CA      4.3219689      1.49       0.96       2.78       3.90       0.11       0.19       0.15
       00017  CB      0.1425789      0.15       0.11       0.06       0.00       0.16       0.02       0.08
       00021  C      19.7719864      1.22       1.07       3.15      16.89       2.66       0.05       0.09
       00022  O       2.7121081      0.05       0.03       1.04       2.09       0.51       0.04       0.03
       00023  N       3.8184184      0.68       1.34       0.78       2.73       1.46       0.06       0.10
       00025  CA      4.1841287      1.24       1.46       1.58       3.77       0.37       0.11       0.16
       00027  CB      0.1262551      0.24       0.24       0.12       0.00       0.39       0.03       0.12
       00031  C      14.0406942      0.99       2.14       0.28      17.62       4.80       0.79       0.10
       00032  O       2.2077275      0.20       0.85       0.00       2.37       2.57       1.12       0.04
       00033  OXT     1.7552869      0.07       0.85       0.00       1.61       1.44       0.59       0.02


Examination of these results reveals the following characteristic features:

*  The improper and torsion componets show relatively large values.
*  Main chain carbonyl carbons show particularly large values for two reasons:
   (1) They are involved in the two improper torsions, C-O-Cα-N and N-H-C-Cα.
   (2) Their atomic volume is small because of the tight packing by the nearby atoms.

