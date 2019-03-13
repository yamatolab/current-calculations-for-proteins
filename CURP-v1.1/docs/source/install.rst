============
Installation
============
The current version of the CURP program supports Mac and Linux systems.

System Requirements
===================

The CURP program runs on the following operating system.

Linux
-----

*  OS Version: CentOS 6.x

Mac
---

*  OS Version: Marvericks or Yosemite with XCode 6.1

Installation
=====================

To build the standard (parallel) version of CURP, it is assumed that 
``Open MPI, gfortran, graphviz, graphviz-dev and Python 2.7.x`` are already installed on your local machine. To install these packages, you can use `yum` for linux, and `MacPorts` or `Homebrew` for Mac.
It is possible to build the serial version of CURP by running ``make`` with the "serial" option.

1. Please download the CURP package from the following URL:

   http://www.comp-biophys.com/yamato-html/curp.html

2. Set the home directory for the CURP package

   If you use the bash interpreter, the setting is done on your machine by editing ``$HOME/.bashrc`` file. For instance, if the CURP directory is $HOME/opt/curp``, append the following line to the .bashrc file. ::

      export CURP_HOME=$HOME/opt/curp

   Then, type in the following command:

      $ source $HOME/.bashrc

3. Change directory to $CURP_HOME, and run ``make`` to build the standard version of CURP::

      $ make

   The CURP scripts chooses the `gfortran` compiler by default. Alternatively, you can build the intel version of the program by running ``make intel`` if you have the Intel Fortran Compiler (ifort) installed ::

      $ make intel

   The serial version of CURP is built by::

      $ make serial

   To clean up the library files and temporary files under ``$CURP_HOME``,
   enter the following command::

      $ make clean

   Note that `hdf5` and `netcdf` source files are downloaded from the Internet by running the `make` command to build CURP. These files remain in the `curp-environ` directory under ``$CURP_HOME`` and are not removed by `make clean`. 
   If you want to delete the curp-environ directory, enter the following command::

      $ make allclean

Running CURP
============

To run CURP, enter the following command:

   $ $CURP_HOME/bin/curp <input_file> > log

The users set the CURP parameters in <input_file> and the messages will be written to ``log``.

For test calculations of CURP, some examples are provided in the ``$CURP_HOME/test`` directory with more examples found in ``$CURP_HOME/tutorial``.

By running ``$CURP_HOME/test/runall.sh``, all of the test examples will be tested.
