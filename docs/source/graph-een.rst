================
curp graph-een
================

.. warning::

    This command is not supported in CURP 3.0.
    Please use the previous version of CURP for this command.

.. code-block:: console

    $ curp graph-een [-h] [-t [TARGETS [TARGETS ...]]]
                     [--forced-output-nodes [FORCE_NODES [FORCE_NODES ...]]]
                     [-p [CLOSE_PAIRS [CLOSE_PAIRS ...]]] [-c CLUSTER_FN]
                     [-n NODE_FN] [--with-one-letter] [-f FIG_FN]
                     [-r THRESHOLD] [-R TPL_THRESHOLD] [--ratio RATIO]
                     [-s GRAPH_SIZE] [--title TITLE] [--direction {LR,TB}]
                     [-I] [--show-negative-values] [--alpha ALPHA]
                     [-lv [LINE_VALUES [LINE_VALUES ...]]]
                     [-lc [LINE_COLORS [LINE_COLORS ...]]]
                     [-lt [LINE_THICKS [LINE_THICKS ...]]]
                     [-lw [LINE_WEIGHTS [LINE_WEIGHTS ...]]]
                     data_fns [data_fns ...]

Show network chart of the energy conductivities.

Positional arguments
--------------------

    `data_fns`
        
        Specify filenames.


Options
-------

    `-h, --help`
    
        Show this help message and exit.

    `-t [TARGETS [TARGETS ...]], --targets [TARGETS [TARGETS ...]]` (defalut: 1-)

        This option allow you to show only target and thier neighbhor nodes. 
        ex), [-t 1-5], [-t 2-9 15 17], [-t 25-].
    
    `--forced-output-nodes [FORCE_NODES [FORCE_NODES ...]]` defalut: [])

        nodes to forcibly show.
        ex), [-t 1-5], [-t 2-9 15 17], [-t 25-].
    
    `-p [CLOSE_PAIRS [CLOSE_PAIRS ...]], --bring-node-pair-close-together [CLOSE_PAIRS [CLOSE_PAIRS ...]]` (defalut: [])

        bring the node pair close together on the EEN graph.
        ex), 1:2, 3:5, 3:15.
    
    `-c CLUSTER_FN, --cluster-filename CLUSTER_FN` (defalut: None)

        Specify a cluster file name.
    
    `-n NODE_FN, --node-style-filename NODE_FN` (defalut: None)

        Specify a file name which node style difinitions are written.
    
    `--with-one-letter` (defalut: False)
        
        use 1 letter representation for amino acid name.
    
    `-f FIG_FN, --output-een-filename FIG_FN`

        Specify a filename to graph EEN.
    
    `-r THRESHOLD, --threshold THRESHOLD` (type: float, defalut: False)
        
        Threshold value to show nodes on chart figure.
    
    `-R TPL_THRESHOLD, --topology-threshold TPL_THRESHOLD` (type: float, defalut: None)
        
        Threshold value to define graph topology. 
        If you don't give this parameter, the value of threshold to show nodes will be used.

    `--ratio RATIO` (type: float, defalut: None)

        Ratio of height/width to chart a figure. 
        If '--graph-size' is set, this option will be overwritten by it.
    
    `-s GRAPH_SIZE, --graph-size GRAPH_SIZE` (defalut: None)

        The graph size in inch unit. 
        ex) -s '3.4,2.5'
    
    `--title TITLE` (defalut: None)
    
        Specify the title of the figure.
    
    `--direction {LR,TB}` (defalut: TB)
    
        Specify the kind of the direction.
    
    `-I, --hide-isolated-nodes`
    
        Hide isolated nodes when appling multiple irEC files.
    
    `--show-negative-values`
    
        show the nodes pair which have negative values.
    
    `--alpha ALPHA` (type: int, defalut: None)
    
        The transparency value of the background color.
    
    `-lv [LINE_VALUES [LINE_VALUES ...]], --line-values [LINE_VALUES [LINE_VALUES ...]]`  (type: float, default: [0.015, 0.008, 0.003])
    
        The threshold values for line attributes. 
        Number of elements in the list must be equal with all line attribute.
    
    `-lc [LINE_COLORS [LINE_COLORS ...]], --line-colors [LINE_COLORS [LINE_COLORS ...]]` (type: str, default: ['red', 'blue', 'green'])
    
        The colors of line. 
        Number of elements in the list must be equal with all line attribute.
    
    `-lt [LINE_THICKS [LINE_THICKS ...]], --line-thicks [LINE_THICKS [LINE_THICKS ...]]` (type: float, default: [4.0, 4.0, 2.5])
    
        The thickness of line. 
        Number of elements in the list must be equal with all line attribute.
    
    `-lw [LINE_WEIGHTS [LINE_WEIGHTS ...]], --line-weights [LINE_WEIGHTS [LINE_WEIGHTS ...]]` (type: float, default: [5.0, 3.0, 1.0])
    
        The weight of line. 
        Number of elements in the list must be equal with all line attribute.
