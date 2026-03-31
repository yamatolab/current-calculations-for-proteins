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


Setting options of `curp compute`
---------------------------------

.. include:: setting.rst

Input and Output files specification
=====================================

NetCDF specification of flux data file
----------------------------------------------

Please visit NetCDF website to use NetCDF format file for more information.

dimensions
~~~~~~~~~~~

   nframe
      The number of frames. This dimension is unlimited.

   npair
      The number of group pairs.

   ncomponent
      The number of components.

      .. note::
         The contents of the ncomponent is determined by the setting of `output.decomp`, `curp.potential` and `curp.nonbonded_method` in the setting file.

         .. include:: components.rst

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

   ncomponent
      The number of components.

      .. note::
         The contents of the ncomponent is determined by the setting of `output.decomp`, `curp.potential` and `curp.nonbonded_method` in the setting file.
         See "NetCDF specification of flux data file" section for more details.


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

Atom group file
~~~~~~~~~~~~~~~

:: 

   [GROUP1]
   (the array of atom indices in "GROUP1")

   [GROUP2]
   (the array of atom indices in "GROUP2")
   ...

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
Each group name should be unique.
Then the range of the constituent atoms are provided. 
You can spacify the range by using `-` symbol.
You can provide multiple data separated by space, empty line, or tab.
Each atom can belong to only one group. If an atom belongs to multiple groups, the error will be raised.

Group pair file
~~~~~~~~~~~~~~~~

::
   
   [GROUP1]
   GROUP1  GROUP2  GROUP3 ...

   [GROUP2]
   GROUP2  GROUP3 ...

When you want to calculate the flux from "GROUP2" to "GROUP1", you can specify the group pair as "[GROUP1]" and "GROUP2" in the group pair file.

For example, when you use the atom group file above, you can specify the group pairs as follows:

::

   [01_ALA_M]
   01_ALA_S  02_ALA_M

   [01_ALA_S]
   03_ALA_S  02_ALA_S  03_ALA_M

   [02_ALA_M]
   03_ALA_S  02_ALA_S  03_ALA_M


Here the first and second lines mean that the two fluxes are calculated and written: flux from 01_ALA_S to 01_ALA_M and the flux from 02_ALA_M to 01_ALA_M.

You cannot specify the same group pair more than once. If you specify the same group pair more than once, the error will be raised.

In the atom group file, the "GROUP1" should be defined earlier than the "GROUP2" when you perform flux calculation from "GROUP2" to "GROUP1". 
In the case of the example above, you cannot specify "01_ALA_M" in the fifth line because "01_ALA_M" is defined earlier than "01_ALA_S" in the atom group file.
If the "GROUP1" is defined later than the "GROUP2", the error will be raised. 

The flux calculation of group pair "GROUP1" and "GROUP1" is available for the intra-group thermal conductivity calculation.
(Note that you cannot calculate the intra-group energy flux because the energy flux should be defined as the flux between different groups.)

Contact
========

Takahisa YAMATO, Dr. Sci. 

Graduate School of Science, Nagoya University,

Furo-cho, Chikusa-ku, Nagoya, 4648602, Japan.

Email: yamato@nagoya-u.jp

https://www.comp-biophys.com/
