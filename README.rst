======================================
CURP: CURrent calculations in Proteins
======================================

**CURP** permits to compute inter-residue flow of energy or heat and atomic stress tensors in a protein, given atomic coordinates and velocity trajectories obtained through molecular dynamics (MD). 
Energy flow data permit to picture an inter-residue Energy Exchange Network as a graph.

Within thermally fluctuating protein molecules under physiological conditions, tightly packed amino acid residues interact with each other through heat and energy exchanges. 
Non-uniform pattern of heat flow in proteins are illustrated and characterized with a theoretical model based on “local heat conductivity” between each residue pair. 
This model demonstrated characteristic features of “hidden dynamic allostery” in PDZ domain [1]_ and allosteric transition in the oxygen sensor domain of FixL [2]_.
Also we applied it to a small protein to understand the features of local thermal transport of protein[3-5]_.

Offical website and tutorial can be found at `<https://curp.jp/>`_.

Installation
============
CURP requires Python3.6 with numpy to work.
You can install python here_, or anaconda there_.

.. _here: https://www.python.org/downloads/release/python-3613/
.. _there: https://www.anaconda.com/download

Install CURP via pip
--------------------

::

    pip install curp

Get CURP from source code 
-------------------------

::

    git clone https://github.com/yamatolab/current-calculations-for-proteins.git
    cd current-calculations-for-proteins
    pip install .

Development
===========
Please read DEVELOP.rst before starting to develop CURP.

References
==========

.. [1] Ishikura, T.; Iwata, Y.; Hatano, T.; Yamato, T. Energy exchange network of inter-residue interactions within a thermally fluctuating protein molecule: A computational study. *J. Comput. Chem.* **2015**, 36:1709-1718

.. [2] Ota, T.; Yamato, T. Energy Exchange Network Model Demonstrates Protein Allosteric Transition: An Application to an Oxygen Sensor Protein. *J. Phys. Chem. B* **2019**, 123:768-775

.. [3] Yamato, T.; Wang, T.; Sugiura, W.; Laprévote, O.; Katagiri, T. Computational study on the thermal conductivity of a protein. *J. Phys. Chem. B* **2022**, 126:3029-3036

.. [4] Wang, T.; Yamato, T.; Sugiura, W; Site-selective heat current analysis of α-helical protein with linear-homopolymer-like model. *J. Phys. Chem. B* **2023**, 158

.. [5] Wang, T.; Yamato, T.; & Sugiura, W. Thermal Energy Transport through Nonbonded Native Contacts in Protein. *J. Phys. Chem. B* **2024**