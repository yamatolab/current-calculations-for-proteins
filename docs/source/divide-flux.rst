==========
curp analyze divide-flux
==========


.. code-block:: console

    curp analyze divide-flux [-h] -o OUTPUT_FN [-d DONOR_LINE]
                             [-a ACCEPTOR_LINE] -t DT [-c COLUMN_LINE]
                             flux_fns [flux_fns ...]

Divide the flux data of all into the flux data for each donor and acceptor.


Positional arguments
--------------------

    `flux_fns`

        Specify filenames in order.
        ex.) flux.dat0001, flux.dat0002, ...


Options
-------

    `-h, --help`
        
        Show this help message and exit.

    `-o OUTPUT_FN, --output-filename OUTPUT_FN`
        
        Specify a filename to output.
    
    `-d DONOR_LINE, --donor DONOR_LINE` (default: None)
    
        Donor names to output to files.
    
    `-a ACCEPTOR_LINE, --acceptor ACCEPTOR_LINE`
    
        Acceptor names to output to files.
    
    `-t DT, --dt DT` (default: 1.0)
    
        Time width[ps] per 1 step.
    
    `-c COLUMN_LINE, --column COLUMN_LINE` (default: None)
    
        Column numbers to output.