==================
Development policy
==================

A few guidelines for developers, like new masters or PhD coming in Yamato's
lab.

Installation
============

You can get curp source code by running ::

    git clone https://gitlab.com/yamato97/current-calculations-for-proteins.git

Then, go in the installed directory and use ::

    pip install .

pip uses *setup.py* file to know what to compile and what commands to install.
*setup.py* works using the packages `distutils`_, `setuptools`_ and
`numpy.distutils`_

.. _distutils: https://docs.python.org/2/library/distutils.html#module-distutils
.. _setuptools: https://setuptools.readthedocs.io/en/latest/
.. _numpy.ditutils: https://docs.scipy.org/doc/numpy/reference/distutils.html

Style
=====
We try, as close as possible, to follow 
`**PEP8** <https://www.python.org/dev/peps/pep-0008/>`. **PEP8** is the 
official writing style in Python, and if you want to become a professional
Python developper you will need to use it, so please go read it before adding
any code on curp. The only lose point is the maximum number of characters on
one line, which can go up to 90ish.
A simple way to use it is to have a plug-in which enforces it and does the
indentations automatically. Please google one for your favorite python text
editor (Atom, vim, emax, they all have one).

**New classes and functions** should **ALWAYS** be written with a docstring. 
Docstrings follow the rules of numpydoc, as described in
`numpydoc docstring guide`_.

You will see that the entirety of CURP doesn't follow these guidelines, as 
these were not enforced before, and lack some substantial commentaries. You are
welcome to add them and "PEP8ify" when reading this old source code if you have
the time, the will and the energy.
.. _numpydoc docstring guide: https://numpydoc.readthedocs.io/en/latest/format.html

API
===
CURP has a shell Application Programming Interface (API). The entire API is
defined in console.py, script/console.py and tools/console.py, and implemented
through setup.py, which is used during "pip install" commands.
This can sometimes provide some information as of what input to give to a
function, if the docstring was not written. If you want to implement a new
command, please pass through this, and if you developed a whole bunch of new
commands, write another console.py and integrate it.

Testing
=======
You wrote some code, good! Now it is time to test it before you add it on
*develop* branch. There are two kinds of test in CURP development:
*unit tests* and *system tests*.

Unit test
---------

Test units are not for the whole system: you just want to test an object or
a function in them. They are pretty important: first, you won't need to reopen
a Python console everytime you *pip install* some code you just implemented and
want to test it; second, if someone writes some code after you and break some 
things, they permit to target what has been broken; finally, they give some
easy and understandable insights to code you wrote for new developers.
You will find test units in the folder they test, in the *test* directory.
With CURP, we use *nose* framework for test units. A good introduction to
*nose* and to testing in Python can be found on the excellent `Python
Testing website`_ from Brian Okken.

Once again, you wiil see a lot of directories not tested, and some tests not
working anymore. Previous students working in Yamato's lab may not have been
aware of good test practices, which is why this happened. Nevertheless, you
are welcome to update any unit test, to write some for old functions, and of
course you have to write unit tests for new ones!

.. _Python Testing website: https://pythontesting.net/framework/nose/nose-introduction/

System test
-----------
System tests use the API to run. Unit test are "low level": they test if a
function or an object works, system tests see if a high level command works,
depending on a lot of objects and functions working together in a real case
scenario.
System tests can be found in *test* directory of Curp GitLab repository. They
should be ran before merging with develop to make sure that nothing is broken.
They should also be updated to consider new curp functionnalities.

These tests are pretty basic, currently if the program executed without error
they will pass, regardless of the result. It would be good to implement some
results comparison: again, feel free to do it if you have the time and will.

Git
===
Git is the versioning system used by CURP on GitLab. It's an important part of
software development: it permits to keep track of changes, hence it permits you
to do some big "ctrl+Z" if you suddenly understand you messed up, and this way
you won't need hundred different folders like "curp-v1.1", "curp-md1",
"curp-md-final". A good introduction to how to use Git in a console can be
found on the website "[Learn Git Branching][https://learngitbranching.js.org/]".

For different projects, you can use different branches, for example if you are
working on a new feature, you make a branch "new-feature" from the branch
"develop", then when you are done with the development and testing of this
new feature, you merge it on develop. 
**New branches** should be made only from develop branch, except for hotfixes
(fixes on a bug on master). Same rule applies for merges.
The develop branch is then merged to master, see `A successful Git branching model`_.

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

.. _A successful Git branching model: https://nvie.com/posts/a-successful-git-branching-model/
.. _How to Write a Git Commit Message: https://chris.beams.io/posts/git-commit/ 


PyPI
----

The `Python Package Index <https://pypi.org/>` (PyPI) is what makes possible
the magic behind pip for installing CURP as a package. CURP has its own
`package page`. To learn more about how to update a package (after pushing
modifications to the master branch and updating the version number), you can
read `Packaging Python projects <https://packaging.python.org/tutorials/packaging-projects/>`
tutorial. Always try first to upload and retrieve the package on test.pypi.org
before updating the true repository.

Building wheels on linux won't work, hence you will either have to use a
continuous integration/deployment server, or simply a Maccintosh.
