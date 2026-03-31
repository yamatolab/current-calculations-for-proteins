============
curp sum-acf
============

.. code-block:: console

   $ sum-acf [-h] -a FILE [-t FILE] [-c COEFFICIENT] ACF_FILE [ACF_FILE ...]

Average auto-correlation function(acf) over the given trajectories.


Positional arguments
---------------------
    
    `ACF_FILE`
        
        ACF file(s) to be averaged.


Options
--------

    `-h, --help`
    
        Show this help message and exit.
    
    `-a FILE, --acf-file FILE` (defalut: None)
    
        The filepath of acf data.

    `-t FILE, --tcs-file FILE` (defalut: None)
    
        The filepath of tc time series data.

    `-c COEFFICIENT, --coefficient COEFFICIENT` (type: float, defalut: 1.0)
    
        Multiply acf by given coefficient.