=============
curp compute
=============

.. code-block:: console

    curp compute [-h] [-v] [-s] [--output-conf-default] [--output-conf-formatted] 
                 [input_]


Launch flux or stress tensor computations given a configuration file.

Positional arguments
--------------------

    `input` (default: run.cfg)
    
        Specify input filenames.

        The contents of the input files are used to determine the type of computation to be performed,
        as well as the parameters and settings for the computation(see :ref:`Reference`).

Options
-------

    `-h, --help`

        Show this help message and exit.

    `-v, --verbose` (default: False)

        Print out informations. 
        This option can be used to get more detailed output during the computation, 
        which can be helpful for debugging or understanding the progress of the calculations.
        
    `-s, --enable_serial` (default: False)
        
        Calculate in serial, don't calculate in parallel.

    `--output-conf-default` (default: False)
        
        Output the config parameters in ini format with default values.

    `--output-conf-formatted` (default: False)
        
        Output the config parameters in rest style with default values.

