CURP: CURrent calculations in Proteins
======================================

CURP permits to compute inter-residue energy or heat flow and atomic stress tensors in a protein, given atomic coordinates and velocity trajectories obtained through molecular dynamics (MD). Energy flow data permit to picture an inter-residue Energy Exchange Network as a graph.

Installation
------------
CURP requires Python2.7 to work. Python3 compatibility has yet to be realized.
You can install python here_, or anaconda there_.

.. _here: https://www.python.org/downloads/release/python-2716/
.. _there: https://www.anaconda.com/distribution/

You can get curp by running ::
    git clone https://gitlab.com/yamato97/current-calculations-for-proteins.git``.

To install it, go in the installed directory and use ::
    pip install .
    
or your favorite python package manager, like ``conda`` or ``pipenv``.

Development
-----------
**New branches** should be made only from development branch, except for hotfixes. Same rule applies for merges. The development branch is then merged to master, see `a successful branching model`_.

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


.. _a successful branching model: https://nvie.com/posts/a-successful-git-branching-model/
.. _How to Write a Git Commit Message: https://chris.beams.io/posts/git-commit/ 
.. _numpydoc docstring guide: https://numpydoc.readthedocs.io/en/latest/format.html
