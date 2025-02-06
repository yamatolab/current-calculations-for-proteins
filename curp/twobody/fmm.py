# -*- coding: utf-8 -*-
import numpy as np


class FMMCalculator:
    
    def __init__(self):
        # from . import lib_fmm
        # self.__mod_fmm = lib_fmm
        pass
        
    def setup(self, info, table, n_crit=10, theta=0.5):
        self.__n_crit = n_crit
        self.__theta = theta
        self.__charge = info['charge']
        self.table = table
        # calculaterにchargeを渡しておく
        
    def initialize(self, particles):
        # t ごとに実行
        # 1. cellの初期化
        # sourceとtargetにわける
        # 
        pass
    
    def run_fmm(self, particles):
        self.initialize(particles)
        
        
        
    

class FMMCellMaker:

    def setup_cell(self):
        """The class for a cell.
    
        Arguments:
            n_crit: maximum number of particles in a leaf cell.
        
        Attributes:
            nleaf (int): number of leaves in the cell
            leaf (array of int): array of leaf index
            nchild (int):  an integer whose last 8 bits is used to keep track 
            of the empty child cells
            child (array of int): array of child index
            parent (int): index of parent cell
            x, y, z (float): coordinates of the cell's center
            r (float): radius of the cell (half of the side length for cubic cell)
            multipole (array of float): multipole array of the cell
        
        """
        
        self.__nleaf = 0                                        # number of leaves
        self.__leaf = np.zeros(self.__n_crit, dtype=np.int)     # array of leaf index
        self.__nchild = 0                                       # binary counter to keep track of empty cells
        self.__child = np.zeros(8, dtype=np.int)                # array of child index
        self.parent = 0                                         # index of parent cell
        self.x = self.y = self.z = 0.                           # center of the cell
        self.r = 0.                                             # radius of the cell
        self.multipole = np.zeros(10, dtype=np.float)           # multipole array

        charge = self.__charge
        phi = self.__phi
        
        
    def set_root_cell(self, particles, rc):
        
        self.__nleaf = 0                                        # number of leaves
        self.__leaf = np.zeros(self.__n_crit, dtype=np.int)     # array of leaf index
        self.__nchild = 0                                       # binary counter to keep track of empty cells
        self.__child = np.zeros(8, dtype=np.int)                # array of child index
        self.parent = 0                                         # index of parent cell
        self.x = rc[0]
        self.y = rc[1]
        self.z = rc[2]
        self.r = max(rc)                                # radius of the cell
        self.multipole = np.zeros(10, dtype=np.float)   # multipole array

    def get_charge(self):
        return self.__charge
        
        
    def build_tree(self, particles, root, n_crit):
    
        """Construct a hierarchical octree to store the particles and return 
        the tree (list) of cells.
        
        Arguments:
            particles: the list of particles.
            root: the root cell.
            n_crit: maximum number of particles in a single cell.
        
        Returns:
            cells: the list of cells.
        
        """
        # set root cell
        cells = [root]       # initialize the cells list
        
        # build tree
        n = len(particles)
        
        for i in range(n):
            
            # traverse from the root down to a leaf cell
            curr = 0
        
            while cells[curr].nleaf >= n_crit:
                cells[curr].nleaf += 1
                octant = (particles[i].x > cells[curr].x) + ((particles[i].y > cells[curr].y) << 1) \
                    + ((particles[i].z > cells[curr].z) << 2)
                
                # if there is no child cell in the particles octant, then create one
                if not cells[curr].nchild & (1 << octant):
                    self.add_child(octant, curr, cells, n_crit)
            
                curr = cells[curr].child[octant]
            
            # allocate the particle in the leaf cell
            cells[curr].leaf[cells[curr].nleaf] = i
            cells[curr].nleaf += 1
            
            # check whether to split or not
            if cells[curr].nleaf >= n_crit:
                self.split_cell(particles, curr, cells, n_crit)
        
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
        cells[c].r = cells[p].r / 2
        cells[c].x = cells[p].x + cells[c].r * ((octant & 1) * 2 - 1)
        cells[c].y = cells[p].y + cells[c].r * ((octant & 2) - 1    )
        cells[c].z = cells[p].z + cells[c].r * ((octant & 4) / 2 - 1)
        
        # establish mutual reference in the cells list
        cells[c].parent = p
        cells[p].child[octant] = c
        cells[p].nchild = (cells[p].nchild | (1 << octant))


    def split_cell(self, particles, p, cells, n_crit):
        
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
            
            octant = (particles[l].x > cells[p].x) + ((particles[l].y > cells[p].y) << 1) \
                + ((particles[l].z > cells[p].z) << 2)   # finds the particle's octant
        
            # if there is not a child cell in the particles octant, then create one
            if not cells[p].nchild & (1 << octant):
                self.add_child(octant, p, cells, n_crit)
            
            # reallocate the particle in the child cell
            c = cells[p].child[octant]
            cells[c].leaf[cells[c].nleaf] = l
            cells[c].nleaf += 1
            
            # check if the child reach n_crit
            if cells[c].nleaf >= n_crit:
                self.split_cell(particles, c, cells, n_crit)


    def get_multipole(self, particles, p, cells, leaves, n_crit):
    
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
        if cells[p].nleaf >= n_crit:
            for c in range(8):
                if cells[p].nchild & (1 << c):
                    self.get_multipole(particles, cells[p].child[c], cells, leaves, n_crit)
        
        # otherwise cell p is a leaf cell
        else:
            # loop in leaf particles, do P2M
            for i in range(cells[p].nleaf):
                l = cells[p].leaf[i]
                dx, dy, dz = cells[p].x-particles[l].x, \
                            cells[p].y-particles[l].y, \
                            cells[p].z-particles[l].z
                
                cells[p].multipole += particles[l].m * \
                                    np.array((1, dx, dy, dz,\
                                                dx**2/2, dy**2/2, dz**2/2,\
                                                dx*dy/2, dy*dz/2, dz*dx/2)) 
            leaves.append(p)


    def M2M(self, p, c, cells):
        
        """Calculate parent cell p's multipole array based on child cell c's 
        multipoles.
        
        Arguments:
            p: parent cell index in cells list.
            c: child cell index in cells list.
            cells: the list of cells.
        """
        
        dx, dy, dz = cells[p].x-cells[c].x, cells[p].y-cells[c].y, cells[p].z-cells[c].z
        
        Dxyz =  np.array((dx, dy, dz))
        Dyzx = np.roll(Dxyz,-1) #It permutes the array (dx,dy,dz) to (dy,dz,dx) 
        
        cells[p].multipole += cells[c].multipole
        
        cells[p].multipole[1:4] += cells[c].multipole[0] * Dxyz
        
        cells[p].multipole[4:7] += cells[c].multipole[1:4] * Dxyz\
                                + 0.5*cells[c].multipole[0] *  Dxyz**2
        
        cells[p].multipole[7:] += 0.5*np.roll(cells[c].multipole[1:4], -1) *  Dxyz \
                                + 0.5*cells[c].multipole[1:4] * Dxyz \
                                + 0.5*cells[c].multipole[0] * Dxyz * Dyzx   

    def upward_sweep(self, cells):
    
        """Traverse from leaves to root, in order to calculate multipoles of all the cells.
        
        Arguments:
            cells: the list of cells.    
        """
    
        for c in range(len(cells)-1, 0, -1):
            p = cells[c].parent
            self.M2M(p, c, cells)

    #----------potential evaluation: particle-particle-----#

    def evaluate(self, particles, p, i, cells, n_crit, theta):
        
        """Evaluate the gravitational potential at a target point i, 
        caused by source particles cell p. If nleaf of cell p is less 
        than n_crit (leaf), use direct summation. Otherwise (non-leaf), loop
        in p's child cells. If child cell c is in far-field of target particle i,
        use multipole expansion. Otherwise (near-field), call the function
        recursively.
        
        Arguments:
            particles: the list of particles
            p: cell index in cells list
            i: target particle index
            cells: the list of cells
            n_crit: maximum number of leaves in a single cell
            theta: tolerance parameter    
        """
    
        # non-leaf cell
        if cells[p].nleaf >= n_crit:
            
            # loop in p's child cells (8 octants)
            for octant in range(8):
                if cells[p].nchild & (1 << octant):
                    c = cells[p].child[octant]
                    r = particles[i].distance(cells[c])
                    
                    # near-field child cell
                    if cells[c].r > theta*r:
                        self.evaluate(particles, c, i, cells, n_crit, theta)
                    
                    # far-field child cell
                    else:
                        dx = particles[i].x - cells[c].x
                        dy = particles[i].y - cells[c].y
                        dz = particles[i].z - cells[c].z
                        r3 = r**3
                        r5 = r3*r**2
                    
                        # calculate the weight for each multipole
                        weight = [1/r, -dx/r3, -dy/r3, -dz/r3, 3*dx**2/r5 - 1/r3, \
                                3*dy**2/r5 - 1/r3, 3*dz**2/r5 - 1/r3, 3*dx*dy/r5, \
                                3*dy*dz/r5, 3*dz*dx/r5]
                        
                        particles[i].phi += np.dot(cells[c].multipole, weight)
                    
        # leaf cell
        else:
            # loop in twig cell's particles
            for l in range(cells[p].nleaf):
                source = particles[cells[p].leaf[l]]
                r = particles[i].distance(source)
                if r != 0:
                    particles[i].phi += source.m / r


    def eval_potential(self, particles, cells, n_crit, theta):
        
        """Evaluate the gravitational potential at all target points 
        
        Arguments:
            particles: the list of particles.
            cells: the list of cells.
            n_crit: maximum number of particles in a single cell.
            theta: tolerance parameter.    
        """
        
        for i in range(len(particles)):
            self.evaluate(particles, 0, i, cells, n_crit, theta)



