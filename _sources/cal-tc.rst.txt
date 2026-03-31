===========
curp cal-tc
===========

.. code-block:: console

    curp cal-tc [-h] -o FILE [-a FILE] [-af FORMAT] [-r FIRST LAST INTER]
                [-s SHIFT] [--sample-number NSAMPLE] [-dt DT] [-c COEFFICIENT]
                [-v]
                FLUX_FILENAME 


Calculate transport coefficients(tc) from energy and heat flux data.


Positional arguments
----------------------

    `FLUX_FILENAME`
    
        The filepath of flux data.


Options
-------

    `-h, --help`

        Show this help message and exit.

    `-o FILE, --tc-file FILE`

        The filename of tc data(output).

    `-a FILE, --acf-file FILE` (default: None)

        The filename of auto correlation function(acf) data.


    `-af FORMAT, --acf-file-format FORMAT {netcdf|ascii}` (default: netcdf)

        The file format of acf data.
        Supported formats are netcdf and ascii.

    `-r FIRST LAST INTER, --frame-range FIRST LAST INTER` (type: int, default: 1, -1, 1)

        The range of flux frame; fist frame, last frame, interval step.

    `-s SHIFT, --average-shift SHIFT` (type: int, default: 1)

        The frame to shift for averaging.

    `--sample-number NSAMPLE` (type: int, default: 0)

        number of sample for one flux data file.
        Default value is 0 that present to make samples as much as possible.

    `-dt DT, --dt DT` (type: float, default: None)

        t of between the neighbour frames. The unit is in ps. 
        Default value is determined by time variable in flux data.

    `-c COEFFICIENT, --coefficient COEFFICIENT` (type: float, default: False)

        Multiply acf by given coefficient.

    `-v, --verbose` (default: False)

        Turn on debug mode.
        
