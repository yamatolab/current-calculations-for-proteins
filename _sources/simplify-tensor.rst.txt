============================
curp analyze simplify-tensor
============================


.. code-block:: console

    curp analyze simplify-tensor [-h] -i FILENAME [-l LABELS] [-s]
                                 [fns [fns ...]]


Show and make figure for stress ratio.


Positional arguments
---------------------

    `fns [fns ...]`
        
        Specify additive filenames.
        ex.) label_data1, label_data2, ...


Options
-------

    `-h, --help`
    
        Show this help message and exit.

    `-i FILENAME, --input-data FILENAME` (defalut: None)
    
        Specify input filename for the stress data.
    
    `-l LABELS, --labels LABELS`
    
        Specify labels of components to analyze.
        ex.) -l "total,bond,angle,..."
    
    `-s, --every-snapshot` (defalut: False)
        
        Specify flag to average the magnitude for every snapshot.

