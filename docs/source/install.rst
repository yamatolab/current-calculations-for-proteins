============
Installation Guide
============
Welcome to the CURP installation guide. This document will help you set up CURP on Mac, 
Linux, and WSL (Windows Subsystem for Linux) systems.

System Requirements
===================

The CURP program has been tested on the following systems:

Linux
-----

*  OS Version: CentOS 6.x, Ubuntu 18.x

Mac
---

*  OS Version: Mavericks or Yosemite with XCode 6.1

CURP runs on Python 3.6 and requires Numpy to be installed.
The system should be working under OpenMPI.

Installation
============

.. note::

   CURP version 1.x requires Numpy version 1.11.2 ~ 1.16.x to build the package before installation.
   Please install Numpy before installation.

Pip Installation
----------------

Run ``pip install Curp`` or any flavor of pip (``pip install --user curp``, ``pipenv install curp``, ...).

If the installation of Curp fails due to the absence of netCDF4, OpenMPI, or other dependencies, 
please install them manually. You can set up the environment using conda, which is the recommended approach.

Conda Installation (Recommended)
--------------------------------

For systems like HPC or WSL where necessary libraries might not be installed by default in the user space, 
you can create a virtual environment with the required dependencies using conda:

.. code-block:: bash

   conda create -n curp-v2 python=3.6 pygraphviz mpi4py netcdf4

Activate the environment and install CURP:

.. code-block:: bash

   conda activate curp-v2
   pip install curp

This method ensures that dependencies, including Fortran and C/C++ compilers, 
are correctly installed and managed within the virtual environment.

Verification
------------

To verify your installation, run the following command:

.. code-block:: bash

   curp -h

This should display the help message of CURP, confirming that the installation was successful.

Latest Development Version
--------------------------
To install the latest development version, run:

.. code-block:: bash

   pip install git+https://github.com/yamatolab/current-calculations-for-proteins.git@develop

Git Installation
----------------

If you want to develop CURP or simply access its source, you can install it through its GitHub repository:

.. code-block:: bash

   git clone https://github.com/yamatolab/current-calculations-for-proteins.git

In the created directory, simply run ``pip install --user .`` or ``python setup.py install --user``. 
If you are an admin and want every user of the machine to be able to use CURP, use ``pip install .`` instead.

Running CURP
============

To run CURP, simply enter ``curp compute <input_file>``.
A few commands are also available from the terminal to analyze CURP results: ``curp cal-tc``, ``curp sum-tc``, ``curp conv-trj`` and ``curp graph-een``.

The users set the CURP parameters in <input_file> (default: run.cfg).

Test calculations of CURP can be run by cloning the GitHub repository, then running ``runall.sh`` in the ``test`` folder.

Troubleshooting
===============

If you encounter any issues during installation, consider the following:

- Ensure all dependencies are installed correctly, especially compilers and libraries.
- Consult the CURP GitHub repository for any known issues or updates.
