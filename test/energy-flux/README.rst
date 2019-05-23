Simple informations of example script

atom
   The target system is ALA-ALA-ALA.
   group_method = united
   flux_grain = atom
   target_atoms = 11-15
   decomp = no

both
   The target system is ALA-ALA-ALA.
   group_method = residue
   flux_grain = both
   target_atoms = 10-15
   decomp = no

residue-decomp
   The target system is ALA-ALA-ALA.
   group_method = residue
   flux_grain = group
   target_atoms = 1-33  
   decomp = yes

residue-pair
   The target system is ALA-ALA-ALA.
   group_method = residue
   flux_grain = group
   target_atoms = 1-33
   group_pair_file = gpair.ndx
   decomp = yes

custom-decomp
   The target system is ALA-ALA-ALA.
   group_method = file
   flux_grain = group
   target_atoms = 11-30
   decomp = yes

water-ala3
   The target system is ALA-ALA-ALA solvated by the water molecules.
   The energy fluxes between residues and solvent are calculated.
   group_method = file
   flux_grain = group
   decomp = no

shell-parallel
   Parallel calculation example with shell script.

b2AR-residue
   The target system is inactive β2 aderenergic receptor with ligand, removed solvation.

   group_method = residue
   flux_grain = group
   group_pair_file = gpair.ndx
   remove_trans = yes
   remove_rotate = no
   decomp = no

pickup
   An example to $CURP_HOME/script/pickup_respairs.py to pickup the residue pairs which are within 5 Å each other.

runall.sh
   The run scripts in all example directories are run.

cleanall.sh
   The output data in all the example directories are deleted.
