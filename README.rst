CURP: CURrent calculations in Proteins
======================================

Installation
------------

CURP requires Python2.7 to work. Python3 compatibility has yet to be realized.
You can install python here_, or anaconda there_.
If you already use python or anaconda but not a compatible version, or even if you start using python, it is advised to rather use pyenv_ for the version gestion.

.. _here: https://www.python.org/downloads/release/python-2716/
.. _there: https://www.anaconda.com/distribution/
.. _pyenv: https://github.com/pyenv/pyenv

Following Python libraries are needed:
 - numpy 1.11.2
 - mpi2py 2.0.0
 - netCDF4 1.2.4
 - nose
 - benchmarker
 - setproctitle
 - epydoc
 - pygraphviz

These requirements are all contained in the requirements.txt. To install them, you can run in the curp directory ``pip install -r requirements.txt`` or if using anaconda ``conda install --yes --file requirements.txt``.

If you want the packages to be installed only for the user, you can run ``pip install --user -r requirements.txt`` 

Installing the packages in a virtual environment:
If using pipenv, you can run ``pipenv install`` in your curp directory. For installing pipenv, ``pip install --user pipenv`` should be sufficient if you already have Python. For more pipenv installation tweaks, please refer to the `official documentation`_.

.. _official documentation: https://pipenv.readthedocs.io/en/latest/install/#installing-pipenv

Alternatively, if you use anaconda, you can run::

conda create -n curp python=2.7 anaconda
source activate curp
conda install -n curp --file requirements.txt

/!\ TODO: Alternatively, if you simply want a totally independent curp, you can run ``make environment``.

Currently there is no setup.py to compile fortran code, you will therefore have to run ``make`` to finally get curp to work.
