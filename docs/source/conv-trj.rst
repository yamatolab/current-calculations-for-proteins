=============
curp conv-trj
=============

.. code-block:: console

   $ curp conv-trj [-h] [-crd | -vel] -i TRJ_FILENAME [TRJ_FILENAME ...] -if FORMAT 
                   [--irange FIRST LAST INTER] [-o TRJ_FILENAME]
                   [-of FORMAT] [--orange FIRST LAST INTER] -pf FORMAT -p TPL_FILENAME 
                   [-nf DIST_FORMAT]
                   {convert-only, dry-run, adjust-vel, mask, dist} ... 


Various scripts to process and analyze trajectories.

Positional arguments
--------------------

    `convert-only`
    
        Convert the trajectory intoother format.

    `dry-run`
    
        Do dry-run mode.

    `adjust-vel`
    
        Adjust time of the velocity trajectory to t from t-dt/2.

    `mask`
    
        Remove solvent from trajectory.

    `dist`
    
        Calculate inter-residue distances.


Common options
--------------

    `-h, --help`
        Show this help message and exit.

    `-crd`
        Specify the format of the trajectory.
        This argument allows you specify multiple formats.

    `-vel`
        Specify the format of the trajectory.
        This argument allows you specify multiple formats.
    
    `-i TRJ_FILENAME [TRJ_FILENAME ...], --input-filenames TRJ_FILENAME [TRJ_FILENAME ...]`
                    
        The trajectory file names.

    `-if FORMAT, --input-formats FORMAT`
        Specify the format of the trajectory.
        This argument allows you specify multiple formats.

    `--irange FIRST LAST INTER`

        The trajectory range to process over all trajectory file.

    `-o TRJ_FILENAME, --output-filename TRJ_FILENAME`
                    
        The trajectory file name.

    `-of FORMAT, --output-format FORMAT`

        The trajectory file format for output.

    `--orange FIRST LAST INTER`
                    
        The trajectory range for output.

    `-pf FORMAT, --topology-format FORMAT`

        Specify the format of the topology file.

    `-p TPL_FILENAME, --topology-file TPL_FILENAME`

        Specify the topology file.

    `-nf DIST_FORMAT, --name-format DIST_FORMAT`

        Specify the format for representing residue identify.


Additional options
----------------

    `convert-only` options
    
        `-h, --help`
        
            Show this help message and exit

    `dry-run` options
    
        No additional options.

    `adjust-vel` options
    
        No additional options.

    `mask` options
            
        No additional options.
    
    `dist` options

        `-m DISTANCE_METHOD, --method DISTANCE_METHOD`
            
            The method used in calculating the inter-residue distances.
        
        `-c CUTOFF_LENGTH, --cutoff-length CUTOFF_LENGTH`
            
            The distance cutoff in angstrom.
