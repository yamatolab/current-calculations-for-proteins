=======================
curp analyze get-ncdata
=======================


.. code-block:: console

   $ curp analyze get-ncdata [-h] -r [FIRST:LAST [FIRST:LAST ...]] [-n DATANAME]
                             [-o PREFIX]
                             ACF_FILE


Get simple text data from file in netcdf format by given arguments.


Positional arguments
---------------------
    
    `ACF_FILE`

        The filepath of auto-correlation function data.


Options
-------

    `-h, --help`
    
        show this help message and exit.

    `-r [FIRST:LAST [FIRST:LAST ...]], --group-ranges [FIRST:LAST [FIRST:LAST ...]]` (default: [])

        The pair range list to want to gain.
    
    `-n DATANAME, --dataname DATANAME` (default: "acf")
    
        The name of the netcdf data you want to gain.
    
    `-o PREFIX, --output-prefix PREFIX` (default: "acf")
    
        The prefix of the files to want to write, that includes directory path.
