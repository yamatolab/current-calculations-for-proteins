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

CURP runs on python2.7 and requires numpy to be installed.
The system should be working under openMPI.

Installation
============

Pip
---

Run ``pip install Curp`` or any flavor of pip (``pip install --user curp``, ``pipenv install curp``, ...) 

Git
---

If you want to develop curp or simply access to its source you can install it through its GitLab repository: ::

git clone https://gitlab.com/yamato97/current-calculations-for-proteins.git

In the created directory, simply run ``pip install --user .`` or ``python setup.py install --user``. If you are admin and want every user of the machine to be able to use curp, use ``pip install .`` instead.

Running CURP
============

To run CURP, simply enter ``curp compute <input_file>``.
A few commands are also available from the terminal to analyze CURP results: ``curp cal-tc``, ``curp sum-tc``, ``curp conv-trj`` and ``curp graph-een``.

The users set the CURP parameters in <input_file> (default: run.cfg).

Test calculations of curp can be run by cloning the GitLab repository, then running ``runall.sh`` in ``test`` folder.
