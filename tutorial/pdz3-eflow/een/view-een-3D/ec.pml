
viewport 520, 520
reset 
hide
bg_color white

# alpha 3 helix
sel a3, (resi 393:399)
create a3_copy, a3
color red, a3_copy
set cartoon_oval_width   , 1.0, a3_copy
set cartoon_oval_length  , 1.5, a3_copy
# set cartoon_transparency, 0.2, a3_copy


sel protein, (resi 306:392 or resi 400:403)
# sel xwater,  (resi 899:906 and name O)
# sel ligand,  resn CAU
# create waters,  (resn WAT and name O)
create ca_protein, name CA

# all
show cartoon
# set cartoon_trace, 1
# spectrum count, rainbow_rev, protein
color gray, protein
set cartoon_transparency, 0.6

# ligand binding pocket
sel lbp, (resi 322,323,324,325,326,327,328,331,339,372,376,379,380)
show surface, lbp
set solvent_radius, 1.0
set surface_color, yellow
set transparency, 0.6

# set ray_shadows, off

turn x, 140
turn z, -80
turn y, 30
zoom ca_protein, 4
