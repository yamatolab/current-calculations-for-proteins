============================
curp analyze pickup-respairs
============================


.. code-block:: console

    $ curp analyze pickup-respairs [-h] -if FORMAT -p PRMTOP_FN -pf PRMTOP_FMT
                                   [-i INTERVAL] [-m {com,nearest,farthest}]
                                   [-c CUTOFF]
                                   [-t [TRIM_RESNAMES [TRIM_RESNAMES ...]]]
                                   [-f FORMAT] [-U]
                                   [-e [FIRST:LAST [FIRST:LAST ...]]] [-b]
                                   trj_fns [trj_fns ...]


Pick up the residue pair table to given file as a ndx format.


Positional arguments
--------------------

    `trj_fns`

        The trajectory file names.


Options
-------

    `-h, --help`
        
        Show this help message and exit.

    `-if FORMAT, --input-format FORMAT`
    
        Specify formats of trajectories.
    
    `-p PRMTOP_FN, --input-prmtop-file PRMTOP_FN`
    
        Specify the prmtop file to be amber format.
    
    `-pf PRMTOP_FMT, --input-prmtop-format PRMTOP_FMT`
    
        Specify the format of the prmtop file.    
    
    `-i INTERVAL, --interval INTERVAL` (type: int, default: 1)
    
        Specify interval step to perform the calculation for trajectory.
    
    `-m {com,nearest,farthest}, --mode {com,nearest,farthest}`
    
        Cutoff method; com, nearest and farthest for residue.
        "com" is the distance between centers of mass of residues.
        "nearest" is the distance between nearest atoms of residues.
        "farthest" is the distance between farthest atoms of residues.
    
    `-c CUTOFF, --cutoff CUTOFF` (default: 5.0)
    
        Specify the cutoff to pick up.
    
    `-t [TRIM_RESNAMES [TRIM_RESNAMES ...]], --trim-resnames [TRIM_RESNAMES [TRIM_RESNAMES ...]]` (default: [])
    
        Residue names for trimming from the target residues.
    
    `-f FORMAT, --output-format FORMAT` (default: {rid:05}_{rname})
    
        Specify the format for representing residue identify.
    
    `-U, --disble-union`
    
        Whether residue pair table is calculated by union or intersection when collecting over trajectory specify

    `-e [FIRST:LAST [FIRST:LAST ...]], --exclude-resids [FIRST:LAST [FIRST:LAST ...]]` (default: [])

        Residues that you want to exclude.
        ex) 5:10 80:150

    `-b, --enable-residues-both-cxcluded`

        Whether to include residue pairs that exclude residues