============
Installation
============
The current version of the CURP program supports Mac and Linux systems.

System Requirements
===================

The CURP program has been tested on the following systems:

Linux
-----

*  OS Version: CentOS 6.x, Ubuntu 18.x

Mac
---

*  OS Version: Marvericks or Yosemite with XCode 6.1

Installation
=====================

1. Please download the CURP package from the following URL:

   http://www.comp-biophys.com/yamato-html/curp.html

Or, if you have access to the repository, run: ::

git clone https://gitlab.com/yamato97/current-calculations-for-proteins.git

2. Compile and install using pip

In the created directory, simply run ``pip install --user .``. If you are admin and want every user of the machine to be able to use curp, use ``pip install .`` instead.

Running CURP
============

To run CURP, simply enter ``curp <input_file>``.
A few commands are also available from the terminal to analyze CURP results: ``cal_tc`` ``conv_trj`` and ``graph_een``.

The users set the CURP parameters in <input_file> (default: run.cfg) and the messages will be written to ``log``.

For test calculations of CURP, some examples are provided in the ``test`` directory with more examples found in ``tutorial``.

By running ``test/runall.sh``, all of the test examples will be tested.
