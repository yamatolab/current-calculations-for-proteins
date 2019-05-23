from __future__ import print_function
import sys


id_names = []
for line in sys.stdin:
    cols = line.split()
    if not cols[0].startswith("ATOM"): continue
    rname = cols[3]
    if cols[3] == 'WAT': continue

    id_names.append((cols[1], cols[2]))


heavies_names = []
atoms_in_heavies = []

index = 0
for id, name in id_names:
    index += 1
    if name[0] in ['C','O','N']:
        heavies_names.append(name)
        atoms = [id]

        for id, name in id_names[index:]:
            if name[0] in ['H']:
                atoms.append(id)
            else:
                break

        atoms_in_heavies.append(atoms)


for name, atoms in zip(heavies_names, atoms_in_heavies):
    print('[{}_{}]'.format(name, atoms[0]))
    print(' '.join(atoms))
    print()



