=========
Reference
=========



.. Test examples
.. =============
.. 
.. There are in tests directories.
.. 

Specification of curp setting file
==================================

This section introduces the description about the setting for the CURP.
The character ``:`` after each of the keyword and the value means value
types allowed. For example:

:Int: integer value.
:Float: floating value.
:Bool: boolean value. yes or no.
:File: specify file path.
:List[type]: list of the value for the type given in ``[`` and ``]``.
:Choice[A|B|C]: ``A``, ``B`` or ``C`` must be chosen.


Setting options
-----------------

.. include:: setting.txt

Input and Output files specification
=====================================

NetCDF specification of energy flux data file
----------------------------------------------

Please visit NetCDF website to use NetCDF format file for more information.

dimensions
~~~~~~~~~~~

nframe
   The number of frames. This dimension is unlimited.

npair
   The number of group pairs.

ncomponent
   The number of components that contains total, bond, angle, torsion, improper torsion, coulomb14, vdw14, coulomb and vdw of energy flux if user turn on the keyword, decomp.

nchar
   The number of character array.

variables
~~~~~~~~~

time(nframe)
   Array of the calculated time.

donors(npair, nchar)
   Array of donor name at i:sup:`th` pair.

acceptors(npair, nchar)
   Array of acceptor name at i:sup:`th` pair.

components(ncomponent, nchar)
   Array of component names

flux(nframe, npair, ncomponent)
   Array flux data.

NetCDF specification of time-correlation data
-----------------------------------------------

Please visit NetCDF website to use NetCDF format file for more information.

dimensions
~~~~~~~~~~~

nframe
   The number of frames.

npair
   The number of group pairs. This dimension is unlimited.

.. note::
   Note that there isn't ncomponent variable in time-correlation data 
   unlike energy flux data file

nchar
   The number of character array.

variables
~~~~~~~~~

time(nframe)
   Array of the calculated time.

donors(npair, nchar)
   Array of donor name at i:sup:`th` pair.

acceptors(npair, nchar)
   Array of acceptor name at i:sup:`th` pair.

acf(npair, nframe)
   Array flux data.

Group file specification
--------------------------

For example, you can separate the main chain and the side chain parts by
using the following specification:

::

   [01_ALA_M]
   1-6 11-12

   [01_ALA_S]
   7-10

   [02_ALA_M]
   13-16 21-22

   [02_ALA_S]
   17-20

   [03_ALA_M]
   23-26 31-33

   [03_ALA_S]
   27-30

The group names are surrounded by `[` and `]`.
Then the range of the constituent atoms are provided. 
You can spacify the range by using `-` symbol.
You can provide 
multiple data saparated by space, empty line, or tab.


How to use the CURP tools
=========================

.. include:: helps.txt

Contact
========

Takahisa YAMATO, Dr. Sci. 

Graduate School of Science, Nagoya University,

Furo-cho, Chikusa-ku, Nagoya, 4648602, Japan.

Email: yamato@nagoya-u.jp

http://www.tb.phys.nagoya-u.ac.jp/~yamato

