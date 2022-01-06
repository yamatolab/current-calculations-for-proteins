======================================
CURP: CURrent calculations in Proteins
======================================

**CURP** permits to compute inter-residue flow of energy or heat and atomic stress tensors in a protein, given atomic coordinates and velocity trajectories obtained through molecular dynamics (MD). Energy flow data permit to picture an inter-residue Energy Exchange Network as a graph.

Within thermally fluctuating protein molecules under physiological conditions, tightly packed amino acid residues interact with each other through heat and energy exchanges. Non-uniform pattern of heat flow in proteins are illustrated and characterized with a theoretical model based on “local heat conductivity” between each residue pair. This model demonstrated characteristic features of “hidden dynamic allostery” in PDZ domain [1]_ and allosteric transition in the oxygen sensor domain of FixL [2]_.

Offical website and tutorial can be found at `<https://curp.jp/>`_.

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
**New branches** should be made only from development branch, except for hotfixes. Same rule applies for merges. The development branch is then merged to master, see `A successful Git branching model`_.

**Commit messages** should follow these rules:

    1. Separate subject from body with a blank line
    2. Limit the subject line to 50 characters
    3. Capitalize the subject line
    4. Do not end the subject line with a period
    5. Use the imperative mood in the subject line
    6. Wrap the body at 72 characters
    7. Use the body to explain what and why vs. how

For example::

    Derezz the master control program

    MCP turned out to be evil and had become intent on world domination.
    This commit throws Tron's disc into MCP (causing its deresolution)
    and turns it back into a chess game.

These rules, example and more explanations can be found on `How to Write a Git Commit Message`_ article from Chris Beams.

**New classes and functions** should **ALWAYS** be written with a docstring. Docstrings follow the rules of numpydoc, as described in `numpydoc docstring guide`_.

**Test units** use nose, although tests haven't been properly configured yet.

References
==========

.. [1] Ishikura, T.; Iwata, Y.; Hatano, T.; Yamato, T. Energy exchange network of inter-residue interactions within a thermally fluctuating protein molecule: A computational study. *J. Comput. Chem.* **2015**, 36:1709-1718
.. [2] Ota, T.; Yamato, T. Energy Exchange Network Model Demonstrates Protein Allosteric Transition: An Application to an Oxygen Sensor Protein. *J. Phys. Chem. B* **2019**, 123:768-775

.. _A successful Git branching model: https://nvie.com/posts/a-successful-git-branching-model/
.. _How to Write a Git Commit Message: https://chris.beams.io/posts/git-commit/ 
.. _numpydoc docstring guide: https://numpydoc.readthedocs.io/en/latest/format.html
