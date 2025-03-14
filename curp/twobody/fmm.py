# -*- coding: utf-8 -*-
import numpy as np
from . import amberbase


class FMMCalculatorBase(amberbase.TwoBodyForceBase):
    
    def __init__(self, charges, gnames_iatoms_pairs, gpair_table):
        self.__n_crit = self.__setting.curp.coulomb_fmm_cell_contains
        self.__theta = self.__setting.curp.coulomb_fmm_theta
        self.__gnames_iatoms_pairs = gnames_iatoms_pairs
        self.__gpair_table = gpair_table
        self.__gnames = []
        self.__mod_fmm.setup(self.__natom, self.__n_crit, self.__theta, charges)
         
    def initialize(self, crd):
        self.__mod_fmm.initialize(crd)
    
####################################################################################################################   

class FMMCellMaker(FMMCalculatorBase):
    
    def make_cells(self, crd):
        
        # build tree
        all_cells = self.build_all_tree(self.__gnames_iatoms_pairs, crd, self.__n_crit)
        
        return all_cells

    
    def setup_cell(self):
        """The class for a cell.
    
        Attributes:
            nleaf (int): number of leaves in the cell
            leaf (array of int): array of leaf index
            nchild (int):  an integer whose last 8 bits is used to keep track 
            of the empty child cells
            child (array of int): array of child index
            parent (int): index of parent cell
            cx, cy, cz (float): coordinates of the cell's center
            r (float): radius of the cell (half of the side length for cubic cell)
            multipole (array of float): multipole array of the cell
        
        """
        
        cell = []
        cell.nleaf = 0                                                  # number of leaves
        cell.leaf = np.zeros(self.__n_crit, dtype=np.int)               # array of leaf index
        cell.nchild = 0                                                 # binary counter to keep track of empty cells
        cell.child = np.zeros(8, dtype=np.int)                          # array of child index
        cell.parent = 0                                                 # index of parent cell
        cell.rc = np.zeros(3)                                           # center of the cell
        cell.cx, cell.cy, cell.cz = cell.rc[0], cell.rc[1], cell.rc[2]  # center of the cell
        cell.r = 0.                                                     # radius of the cell
        cell.multipole = np.zeros((10), dtype=np.float)                 # multipole array

        return cell
    
    def build_all_tree(self, group_atoms, crd, n_crit):
        all_cells = []
        for gname, atoms in group_atoms:
            root_cell = self.setup_cell(atoms)
            
            all_cells.append(gname, self._build_tree(atoms, crd, root_cell, n_crit))
            self.__gnames.append(gname)
            
        return all_cells
    
        
    def _build_tree(self, atoms, crd, root, n_crit):
    
        """Construct a hierarchical octree to store the particles and return 
        the tree (list) of cells.
        
        Arguments:
            s: the list of particles.
            root: the root cell.
            n_crit: maximum number of particles in a single cell.
        
        Returns:
            cells: the list of cells.
        
        """
        # set root cell
        cells = root       # initialize the cells list
        root.cx, root.cy, root.cx, root.r = self.__mod_fmm.calculate_rc(atoms)
        
        # build tree
        n = len(atoms)
        
        for i in range(n):
            
            # traverse from the root down to a leaf cell
            curr = 0
            np = atoms[i] - 1  # particle`s crd index
            
            while cells[curr].nleaf >= n_crit:
                cells[curr].nleaf += 1
                octant = (crd[np,0] > cells[curr].cx) + ((crd[np,1] > cells[curr].cy) << 1) \
                    + ((crd[np,2] > cells[curr].cz) << 2)
                
                # if there is no child cell in the particles octant, then create one
                if not cells[curr].nchild & (1 << octant):
                    self.add_child(octant, curr, cells, n_crit)
            
                curr = cells[curr].child[octant]
            
            # allocate the particle in the leaf cell
            cells[curr].leaf[cells[curr].nleaf] = i
            cells[curr].nleaf += 1
            
            # check whether to split or not
            if cells[curr].nleaf >= n_crit:
                self.split_cell(crd, np, curr, cells, n_crit)
        
        return cells


    def add_child(self, octant, p, cells, n_crit):
        
        """Add a cell to the end of cells list as a child of p, initialize the
        center and radius of the child cell c, and establish mutual reference
        between child c and parent p.
        
        Arguments:
            octant: reference to one of the eight divisions in three dimensions.
            p: parent cell index in cells list.
            cells: the list of cells.
            n_crit: maximum number of particles in a leaf cell.
        """
        
        # create a new cell instance
        cells.append(self.setup_cell(n_crit))
        
        # the last element of the cells list is the new child c
        c = len(cells) - 1
        
        # geometry relationship between parent and child
        cells[c].r  = cells[p].r / 2
        cells[c].cx = cells[p].cx + cells[c].r * ((octant & 1) * 2 - 1)
        cells[c].cy = cells[p].cy + cells[c].r * ((octant & 2) - 1    )
        cells[c].cz = cells[p].cz + cells[c].r * ((octant & 4) / 2 - 1)
        
        # establish mutual reference in the cells list
        cells[c].parent = p
        cells[p].child[octant] = c
        cells[p].nchild = (cells[p].nchild | (1 << octant))


    def split_cell(self, crd, np, p, cells, n_crit):
        
        """Loop in parent p's leafs and reallocate the particles to subcells. 
        If a subcell has not been created in that octant, create one using add_child. 
        If the subcell c's leaf number exceeds n_crit, split the subcell c recursively.
        
        Arguments: 
            particles: the list of particles.
            p: parent cell index in cells list.
            cells: the list of cells.
            n_crit: maximum number of particles in a leaf cell.
        """
        
        # loop in the particles stored in the parent cell that you want to split
        for l in cells[p].leaf:
            
            octant = (crd[l,0] > cells[p].cx) + ((crd[l,1] > cells[p].cy) << 1) \
                    + ((crd[l,2] > cells[p].cz) << 2)   # finds the particle's octant
        
            # if there is not a child cell in the particles octant, then create one
            if not cells[p].nchild & (1 << octant):
                self.add_child(octant, p, cells, n_crit)
            
            # reallocate the particle in the child cell
            c = cells[p].child[octant]
            cells[c].leaf[cells[c].nleaf] = l
            cells[c].nleaf += 1
            
            # check if the child reach n_crit
            if cells[c].nleaf >= n_crit:
                self.split_cell(crd, c, cells, n_crit)


###################################### Calculator ##########################################################################

#TODO
class FMMCellCalculator(FMMCalculatorBase):

    def __init__(self):
        self.__leaves = []
        
    def cal_fmm(self, all_cells, crd):
        
        # get multipole arrays
        multipole = [self.get_multipole(crd, 0, all_cells[gname], self.__leaves) for gname, atoms in self.__gnames_iatoms_pairs]
        
        # upward sweep
        m2m = [self.cal_M2M(all_cells[i]) for i in all_cells]

        # evaluate potential
        coulomb_fmm = [self.eval_potential(group_i, groups, all_cells) for group_i, groups in self.__gpair_table]
        
        return coulomb_fmm
        

    def get_multipole(self, crd, p, cells, leaves):
    
        """Calculate multipole arrays for all leaf cells under cell p. If leaf
        number of cell p is equal or bigger than n_crit (non-leaf), traverse down
        recursively. Otherwise (leaf), calculate the multipole arrays for leaf cell p.
        
        Arguments:
            p: current cell's index.
            cells: the list of cells.
            leaves: the array of all leaf cells.
            n_crit: maximum number of particles in a leaf cell.     
        """
    
        # if the current cell p is not a leaf cell, then recursively traverse down
        if cells[p].nleaf >= self.__n_crit:
            for c in range(8):
                if cells[p].nchild & (1 << c):
                    self.get_multipole(crd, cells[p].child[c], cells, leaves)
        
        # otherwise cell p is a leaf cell
        else:
            # loop in leaf particles, do P2M
            cells[p].multipole += self.__mod_fmm.cal_multipole(cells[p].multipole, cells[p].rc, cells[p].nleaf) 
            leaves.append(p)


    def cal_M2M(self, cells):
        return self.__mod_fmm.M2M(cells)

    def eval_potential(self, group_i, groups, cells):
        
        """Evaluate the gravitational potential at all target points 
        
        Arguments:
            particles: the list of particles.
            cells: the list of cells.
            n_crit: maximum number of particles in a single cell.
            theta: tolerance parameter.    
        """

        atomwise, fmm = self.__mod_fmm.evaluate(group_i, groups, cells)            
        return dict(atomwise=atomwise, fmm=fmm)


####################################################################################################################

def check_setting(setting):
    """Check setting parameters for FMM calculation."""

    if setting.curp.method != 'energy-flux' or 'heat-flux':
        raise ValueError('The method should be "energy-flux" or "heat-flux"')
    else:
        if setting.curp.flux_grain != 'group':
            raise ValueError('The flux grain should be "group"')
        else:
            return True

