import os, sys

class Atom:
    terminal = False

    def __str__(self):
        return '<Atom, id={}, name={}, rid={}, rname={}>'.format(
                self.id, self.name, self.rid, self.rname)

class System:

    def __init__(self, pdb_fpfd):

        atoms = self._gen_pdbatom(pdb_fpfd)
        self.__atoms = list(atoms)
        self.__natom = len(self.__atoms)
        self.__rid_to_atoms = self._make_rid_to_atoms()

        self.__cur_iatm = 0

    def __next__(self):
        self.__cur_iatm += 1
        if self.__cur_iatm <= self.__natom:
            return self.__atoms[self.__cur_iatm-1]
        else:
            raise StopIteration
    next = __next__
    
    def __iter__(self):
        return self

    def atom(self, iatm):
        return self.__atoms[iatm-1]

    def to_string(self):
        return "\n".join(self.gen_pdbline())

    def residue(self, resid):
        return self.__rid_to_atoms.get(resid)

    def gen_pdbline(self):
        """Generate the pdb string from atom list."""

        line_template = (
            'ATOM  {id:>5} {name:<4} {rname:>3} {cname:>1}{rid:>5}'
            '   {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfact:6.2f}'
            '          {elem:<2}  '
        )

        for iatm_1, atom in enumerate(self.__atoms):
            x, y, z = atom.pos
            yield line_template.format(
                id=iatm_1+1, name=atom.name, rname=atom.rname, cname=atom.chain,
                rid=atom.rid, x=x, y=y, z=z,
                occ=atom.occupancy, bfact=atom.bfactor, elem=atom.element)

            if atom.terminal:
                yield 'TER'

    def _gen_pdbatom(self, pdb_fpfd, vervose=False):
        """Generate atoms from pdb file."""

        rid_shift = 0

        if type(pdb_fpfd) == file:
            pdb_file = pdb_fpfd

        else:
            pdb_fpfd = os.path.expanduser(pdb_fpfd)
            if pdb_fpfd.endswith('.gz'):
                import gzip
                pdb_file = gzip.open(pdb_fpfd, 'r')
            else:
                pdb_file = open(pdb_fpfd, 'r')

        atom_prev = None

        # get first atom
        for line in pdb_file:
            if line[0:6].strip() in ['ATOM', 'HETATM']:
                atom_prev = parse_pdbline(line)
                natm = 1
                break
            else:
                continue

        # get all atoms
        for line in pdb_file:

            if line[0:6].strip() == 'TER':
                atom_prev.terminal = True

            if line[0:6].strip() not in ['ATOM', 'HETATM']: continue

            atom = parse_pdbline(line)
            natm += 1

            yield atom_prev

            # correction of residue id for pdb file obtained by Amber tools.
            if atom_prev.rid == 9999 and atom.rid == 0:
                rid_shift += 10000

            atom.rid += rid_shift
            
            atom_prev = atom

        else:
            atom_prev.terminal = True
            yield atom_prev

        if vervose:
            print('REMARK the number of atoms that was read: {}'.format(natm))

        if pdb_file.name != '<stdin>':
            pdb_file.close()
        
    def _make_rid_to_atoms(self):

        rid_to_atoms = {}
        prev_rid = 1
        atoms = []

        for atom in self.__atoms:

            if atom.rid != prev_rid:
                rid_to_atoms[prev_rid] = atoms
                atoms = []

            atoms.append(atom)
            prev_rid = atom.rid

        else:
            rid_to_atoms[prev_rid] = atoms

        return rid_to_atoms

def parse_pdbline(line):
    atom = Atom()

    try:

        atom.id = int(line[6:11])
        atom.name = line[12:16].strip()
        atom.rname = line[17:20].strip()
        atom.chain = line[21:22].strip()
        atom.rid = int(line[22:28])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        atom.pos = x, y, z

        occ = line[54:60].strip()
        atom.occupancy = float(occ) if occ else 1.00

        bfact =  line[60:66].strip()
        atom.bfactor = float(bfact) if bfact else 0.00
    except:
        print(line)
        raise

    # get elements
    if len(atom.name) == 1:
        atom.element = atom.name[0]
    else:
        atom.element = atom.name[1] if atom.name[0].isdigit() else atom.name[0]

    return atom

def gen_pdbatom(pdb_fn, vervose=False):
    """Generate atoms from pdb file."""

    if pdb_fn.endswith('.gz'):
        import gzip
        pdb_file = gzip.open(pdb_fn, 'rb')
    else:
        pdb_file = open(pdb_fn, 'r')

    atom_prev = None

    # get first atom
    for line in pdb_file:
        if line[0:6].strip() in ['ATOM', 'HETATM']:
            atom_prev = parse_pdbline(line)
            natm = 1
            break
        else:
            continue

    # get all atoms
    for line in pdb_file:

        if line[0:6].strip() == 'TER':
            atom_prev.terminal = True

        if line[0:6].strip() not in ['ATOM', 'HETATM']: continue

        atom = parse_pdbline(line)
        natm += 1

        yield atom_prev
        
        atom_prev = atom

    else:
        atom_prev.terminal = True
        yield atom_prev

    if vervose:
        print('REMARK the number of atoms that was read: {}'.format(natm))

    pdb_file.close()

def getPDB(pdb_fn):
    """ get molecules data from pdb file
    """

    pdb_file = open(pdb_fn, 'r')

    system = []
    for line in pdb_file:
        if line[0:6].strip() != "ATOM": continue
        atom = {}
        atom["id"] = int(line[6:11])
        atom["name"] = line[12:16].strip()
        atom["rname"] = line[17:20].strip()
        atom["rid"] = int(line[22:26])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        atom["xyz"] = [x, y, z]

        system.append(atom)

    pdb_file.close()
    return system


def gen_pdbline(atoms):
    """Generate the pdb string from atom list."""

    line_template = (
            '{record:<6}{id:>5} {name:<4} {rname:>3} {cname:>1}{rid:>5}'
            '   {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfact:6.2f}'
            '          {elem:<2}  '
    )

    for iatm, atom in enumerate(atoms,1):
        x, y, z = atom.pos
        aid = atom.id if hasattr(atom, 'id') else iatm
        yield line_template.format(record='ATOM',
            id=aid, name=atom.name, rname=atom.rname, cname=atom.chain,
            rid=atom.rid%10000, x=x, y=y, z=z,
            occ=atom.occupancy, bfact=atom.bfactor, elem=atom.element)

        if atom.terminal:
            yield 'TER'

def gen_residues(atoms):

    try:
        atom = atoms.next()
    except AttributeError:
        atom  = atoms[0]
        atoms = atoms[1:]

    res_atoms = [atom]
    res_id = atom.rid

    for atom in atoms:
        if res_id != atom.rid:
            yield res_atoms
            res_atoms = []

        res_atoms.append(atom)
        res_id = atom.rid

    yield res_atoms


if __name__ == '__main__':

    fn = '~/Desktop/complex.pdb'
    pdb_fn = os.path.expanduser(fn)

    atoms = gen_pdbatom(pdb_fn)
    for line in gen_pdbline(atoms):
        print(line)

    pdb_fn = './testdata/system.pdb'
    system = System(pdb_fn)
    # print(system.to_string())
    
    # iatm = 5
    # print(system.atom(iatm))
    for atom in system.residue(5):
        print(atom)

    for atom in system.residue(6):
        print(atom)

    for atom in system.residue(1500):
        print(atom)

