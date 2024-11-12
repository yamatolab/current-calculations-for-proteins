======================================
CURP: CURrent calculations in Proteins
======================================

.. image:: https://pepy.tech/badge/curp
    :target: https://pepy.tech/project/curp

**CURP** permits to compute inter-residue flow of energy or heat and atomic stress tensors in a protein, given atomic coordinates and velocity trajectories obtained through molecular dynamics (MD). Energy flow data permit to picture an inter-residue Energy Exchange Network as a graph.

Within thermally fluctuating protein molecules under physiological conditions, tightly packed amino acid residues interact with each other through heat and energy exchanges. Non-uniform pattern of heat flow in proteins are illustrated and characterized with a theoretical model based on “local heat conductivity” between each residue pair. This model demonstrated characteristic features of “hidden dynamic allostery” in PDZ domain [1]_ and allosteric transition in the oxygen sensor domain of FixL [2]_.

Offical website and tutorial can be found at `<http://www.comp-biophys.com/yamato-lab/curp.html>`_.

Installation
============
CURP requires Python2.7 with numpy to work. Python3 compatibility has yet to be realized.
You can install python here_, or anaconda there_.

.. _here: https://www.python.org/downloads/release/python-2716/
.. _there: https://www.anaconda.com/distribution/

To install curp, run ::

    pip install Curp

or your favorite python package manager, like ``conda`` or ``pipenv``.
You can get curp source code by running ::

    git clone https://gitlab.com/yamato97/current-calculations-for-proteins.git

Then, go in the installed directory and use ::

    pip install .

Development
===========
Please read DEVELOP.rst before starting to develop CURP.

References
==========

.. [1] Ishikura, T.; Iwata, Y.; Hatano, T.; Yamato, T. Energy exchange network of inter-residue interactions within a thermally fluctuating protein molecule: A computational study. *J. Comput. Chem.* **2015**, 36:1709-1718
    [`CrossRef1 <https://doi.org/10.1002/jcc.23989>`_]

.. [2] Ota, T.; Yamato, T. Energy Exchange Network Model Demonstrates Protein Allosteric Transition: An Application to an Oxygen Sensor Protein. *J. Phys. Chem. B* **2019**, 123:768-775
    [`CrossRef2 <https://doi.org/10.1021/acs.jpcb.8b10489>`_]
