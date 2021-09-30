==========
About CURP
==========

Underlying Concepts of CURP
============================

CURrent calculation for Proteins (CURP) permits to compute inter-residue energy or heat conductivity and atomic stress tensors in a protein, given atomic coordinates and velocity trajectories obtained through molecular dynamics (MD).

 The physical properties of condensed matter systems can be illustrated at the microscopic scale by evaluating the continuum variables in atomistic simulations. Among condensed matter systems, native protein molecules are particularly interesting in that they perform biological functions via an energy relaxation process from the non-equilibrium state to the equilibrium state in response to external perturbations such as light illumination and ligand binding. We can characterize native proteins in terms of continuum variables and provide biophysical grounds for the molecular mechanisms of protein functions. The values of the continuum variables vary from site to site in the protein molecule of interest, reflecting its inhomogeneity and anisotropy. Therefore, we can expect to find functionally important sites by analyzing the 3D field of the continuum variables in a protein molecule. For instance, we analyzed electron-tunneling currents in proteins and discovered a key residue for the electron transfer reaction in DNA photolyases[1]. Other examples are the analysis of the energy-transfer pathway in a photoactive yellow protein[2,3], and the identification of the “epicenter” of the protein quake in a photoactive yellow protein via stress tensor analysis[4]. Note that the first quantity is related to the electron flow, the second is related to the energy flow, and the third is related to the linear momentum flow.

 The CURrent calculation for Proteins (CURP) program is written in Python and FORTRAN; it reads (1) the parameters of the force-field functions and the molecular topology data and (2) the atomic coordinates and velocities from the molecular dynamics trajectory. These are then used to calculate the atomic stress tensors. The AMBER format is supported for the current version[5]. Note that it is possible to convert Gromacs topology and trajectory files to Amber format with the aid of external programs such as ParmEd and VMD (see CURP top page).

References:

    [1] Y Miyazawa, H Nishioka, K Yura, T Yamato, *Biophys. J.* 94 (2008) 2194-203.

    [2] T Ishikura, T Yamato, *Chem. Phys. Lett.* 432 (2006) 533-37.

    [3] T Yamato, in: D.M. Leitner, J.E. Straub (Eds.), *Proteins : energy, heat and signal flow*, Taylor and Francis, New York, 2009, pp.129-47.

    [4] K Koike, K Kawaguchi, T Yamato, *Phys Chem Chem Phys* 10 (2008) 1400-5.

    [5] T Ishikura, T Hatano, T Yamato, *Chem. Phys. Lett.* 539 (2012) 144-50.

Flow diagram of CURP calculations
=================================

 The CURP program is designed for the analysis of (1) energy flow during vibrational energy relaxation and (2) atomic stress tensors or the flow of linear momenta in the thermally fluctuating protein media. The flowchart of the CURP calculation procedure is shown below.

.. image:: flowchart.png
   :scale: 80%




License
=======

.. literalinclude:: ../../LICENSE.txt
   :tab-width: 5

Citation
========

Please cite the following references when you use the CURP program.

(1) The main citation of the CURP program

T.Ishikura, Y.Iwata, T.Hatano and T.Yamato, "Energy exchange network of inter-residue interactions within a thermally fluctuating protein molecule: A computational study" *J. Comput. Chem. 36:1709-1718 (2015).

(2) Atomic stress tensor analysis

T.Ishikura, T.Hatano and T.Yamato,  "Atomic  stress  tensor analysis of proteins", *Chem. Phys. Lett.* 539:144-150 (2012).

(3) Energy flow analysis of proteins

T.Yamato, "Energy flow pathways in photoactive yellow proteins", in *Proteins: energy, heat, and signal flow*, eds. D.M. Leitner & J.E. Straub, pp. 129-147, (2009), Taylor and Francis.

T.Ishikura and T.Yamato, "Energy transfer pathways relevant for long-range intramolecular signaling of photosensory protein revealed by microscopic energy conductivity", *Chem. Phys. Lett.* 432:533-537, (2006).

D.Leitner and T. Yamato, "Locating energy transport networks in proteins" (to be submitted)
