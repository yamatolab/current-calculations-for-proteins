#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <queue>
#include <functional>
#include <chrono>
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
namespace py = pybind11;
using namespace Eigen;

class cal_fmm{
        
// public variables used for calculation
public:

    int   natom;                    // number of atoms
    int   ngrp;                     // number of group
    int   n_crit;                   // n_crit (max number of atoms included in the smallest cell)
    float theta;                    // theta  (parameter to decide whether to use multipole expansion or not)

    std::vector<double> charges;    // charges of atoms
    MatrixXd t_crd;                 // coordinates of atoms
    MatrixXd t_vel;                 // velocities of atoms

    std::vector<std::pair<int, int>> bonded_pairs;                              // list of bonded pairs
    std::vector<std::pair<std::string, std::vector<std::string>>> gpair_table;  // table of target group pairs
    std::vector<std::pair<std::string, std::vector<int>>> gname_iatoms_pairs;   // list of group name and its atom indices
    VectorXi iatom_to_igroup;                                                   // mapping from atom index to group index

    // vars for vdw
    std::vector<int> atom_types;
    MatrixXd c6s;
    MatrixXd c12s;

    // vectors for heat flux
    std::vector<std::vector<Vector3d>> hflux_ij;
    std::vector<std::vector<Vector3d>> hflux_ij_near;
    std::vector<std::vector<Vector3d>> hflux_ij_far;

    // Matrices for energy flux
    MatrixXd eflux_ij;
    MatrixXd eflux_ij_near;
    MatrixXd eflux_ij_far;

    // flags for flux type
    bool flag_heat;
    bool flag_energy;

    // function for calculating flux
    std::function<void(Vector3d, Vector3d, Vector3d, int, int)> cal_flux_near;
    std::function<void(Vector3d, Vector3d, Vector3d, int, int)> cal_flux_near_bonded;
    std::function<void(Vector3d, Vector3d, Vector3d, Matrix3d, int, int)> cal_flux_far;

    // near and far pairs for each source atom
    std::vector<std::pair<int, std::vector<int>>> pairs_near;   // [[1, [atom1, atom2, ...]], [2, [atom1, atom2, ...]], ...]
    std::vector<std::pair<int, std::vector<int>>> pairs_far;    // [[1, [cell1, cell2, ...]], [2, [cell1, cell2, ...]], ...]

    // constructor
    cal_fmm():
        natom(0),
        ngrp(0),
        n_crit(0),
        theta(0.0),
        charges(std::vector<double>()),
        t_crd(MatrixXd::Zero(natom, 3)),
        t_vel(MatrixXd::Zero(natom, 3)),
        bonded_pairs(std::vector<std::pair<int, int>>()),
        gpair_table(std::vector<std::pair<std::string, std::vector<std::string>>>()),
        gname_iatoms_pairs(std::vector<std::pair<std::string, std::vector<int>>>()),
        iatom_to_igroup(VectorXi::Zero(0)),
        atom_types(std::vector<int>()),
        c6s(MatrixXd::Zero(0, 0)),
        c12s(MatrixXd::Zero(0, 0)),
        hflux_ij(std::vector<std::vector<Vector3d>>()),
        hflux_ij_near(std::vector<std::vector<Vector3d>>()),
        hflux_ij_far(std::vector<std::vector<Vector3d>>()),
        eflux_ij(MatrixXd::Zero(0, 0)),
        eflux_ij_near(MatrixXd::Zero(0, 0)),
        eflux_ij_far(MatrixXd::Zero(0, 0)),
        flag_heat(false),
        flag_energy(false),
        pairs_near(std::vector<std::pair<int, std::vector<int>>>()),
        pairs_far(std::vector<std::pair<int, std::vector<int>>>())
    {};

    // for debugging
    int count_near = 0;     // counter of atom-atom pairs
    int count_far  = 0;     // counter of atom-cell pairs

    bool output_to_err_file = false; // if true, output process time for each function to error file

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Structure of a cell
    struct Cell {
        int               nleaf;        // number of atoms(leaf) in the cell
        std::vector<int>  leaf;         // index of atoms in the cell
        int               nchild;       // number of child cells
        std::vector<int>  child;        // index of 8 child cells
        int               parent;       // index of parent cell
        Vector3d          rc;           // center of the cell
        double            r;            // radius of the cell
        VectorXd          multipole;    // 10 multipoles
        MatrixXd          multipole_j;  // 3x10 multipoles_j

        Cell():
            nleaf(0),
            leaf(std::vector<int>(0)),
            nchild(0),
            child(std::vector<int>(8, 0)),
            parent(0),
            rc(Vector3d::Zero()),
            r(0.0),
            multipole(VectorXd::Zero(10)),
            multipole_j(MatrixXd::Zero(3, 10))
        {}
    };

    // Structure of all cells for each group
    struct All_cells {
        std::string group;
        std::vector<Cell> cells;

        All_cells():
            group(""),
            cells()
        {}
    };
    std::vector<All_cells> all_cells;

    // Structure for calculating the radius and center of the root cell
    struct Root_r{
        Vector3d rc;
        double r;

        Root_r():
            rc(Vector3d::Zero()),
            r(0.0)
        {}
    };
    Root_r root_r;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // setup parameters
    void setup(const int& input_natom, const int& input_n_crit, const float& input_theta, const std::vector<double>& input_charges, \
        const std::vector<std::pair<int, int>>& input_bonded_pairs, \
        const std::vector<std::pair<std::string, std::vector<int>>>& input_gname_iatoms_pairs, \
        const std::vector<std::pair<std::string, std::vector<std::string>>>& input_gpair_table, \
        const VectorXi& input_iatom_to_igroup,
        const std::vector<int>& input_atom_types,
        const MatrixXd& input_c6s,
        const MatrixXd& input_c12s){

        int t0 = time_now();

        // set parameters
        natom   = input_natom;
        n_crit  = input_n_crit;
        theta   = input_theta;
        charges = input_charges;

        bonded_pairs       = input_bonded_pairs;
        gname_iatoms_pairs = input_gname_iatoms_pairs;
        gpair_table        = input_gpair_table;
        iatom_to_igroup    = input_iatom_to_igroup;
        
        atom_types = input_atom_types;
        c6s        = input_c6s;
        c12s       = input_c12s;

        ngrp = iatom_to_igroup.maxCoeff();

        int t1 = time_now();
        if (output_to_err_file == true){
            std::cerr << "setup time: " << t1 - t0 << " seconds" << std::endl;
        }
    };

    // get flux type
    void set_flux(const std::string& flux_type){

        if (flux_type == "heat"){
            flag_heat = true;

            // initialize heat flux variables
            hflux_ij.resize(ngrp, std::vector<Vector3d>(ngrp, Vector3d::Zero()));
            hflux_ij_near.resize(ngrp, std::vector<Vector3d>(ngrp, Vector3d::Zero()));
            hflux_ij_far.resize(ngrp, std::vector<Vector3d>(ngrp, Vector3d::Zero()));

            // determine flux calculation function for heat flux
            cal_flux_far  = [=](const Vector3d fiJ, const Vector3d vi, const Vector3d r, const Matrix3d multi_j, const int igrp, const int Jgrp) {
                cal_hflux_far(fiJ, vi, r, multi_j, igrp, Jgrp);
            };
            cal_flux_near = [=](const Vector3d fij, const Vector3d vi, const Vector3d rij, const int idx_source, const int idx_target) {
                cal_hflux_near(fij, vi, rij, idx_source, idx_target);
            };
            cal_flux_near_bonded = [=](const Vector3d fij, const Vector3d vij, const Vector3d rij, const int igrp, const int jgrp) {
                cal_hflux_near_bonded(fij, vij, rij, igrp, jgrp);
            };

        }
        else if (flux_type == "energy"){
            flag_energy = true;

            // initialize energy flux variables
            eflux_ij      = MatrixXd::Zero(ngrp, ngrp);
            eflux_ij_near = MatrixXd::Zero(ngrp, ngrp);
            eflux_ij_far  = MatrixXd::Zero(ngrp, ngrp);

            // determine flux calculation function for energy flux
            cal_flux_far  = [=](const Vector3d fiJ, const Vector3d vi, const Vector3d r, const Matrix3d multi_j, const int igrp, const int Jgrp) {
                cal_eflux_far(fiJ, vi, r, multi_j, igrp, Jgrp);
            };
            cal_flux_near = [=](const Vector3d fij, const Vector3d vi, const Vector3d rij, const int idx_source, const int idx_target) {
                cal_eflux_near(fij, vi, rij, idx_source, idx_target);
            };
            cal_flux_near_bonded = [=](const Vector3d fij, const Vector3d vij, const Vector3d rij, const int igrp, const int jgrp) {
                cal_eflux_near_bonded(fij, vij, rij, igrp, jgrp);
            };
        }
        else{
            std::cerr << "Invalid flux type" << std::endl;
            std::exit(1);
        }
    };

    // read trajectory, initialize variables for flux calculation
    void initialize(const MatrixXd& crd, const MatrixXd& vel){
        int t0 = time_now();

        // initialize variables
        t_crd = MatrixXd::Zero(natom, 3);
        t_vel = MatrixXd::Zero(natom, 3);
        t_crd = crd;
        t_vel = vel;

        count_near = 0;
        count_far  = 0;

        // initialize flux variables
        if (flag_heat == true){
            hflux_ij.clear();
            hflux_ij_near.clear();
            hflux_ij_far.clear();
            hflux_ij.resize(ngrp, std::vector<Vector3d>(ngrp, Vector3d::Zero()));
            hflux_ij_near.resize(ngrp, std::vector<Vector3d>(ngrp, Vector3d::Zero()));
            hflux_ij_far.resize(ngrp, std::vector<Vector3d>(ngrp, Vector3d::Zero()));
        }
        else if (flag_energy == true){
            eflux_ij      = MatrixXd::Zero(ngrp, ngrp);
            eflux_ij_near = MatrixXd::Zero(ngrp, ngrp);
            eflux_ij_far  = MatrixXd::Zero(ngrp, ngrp);
        }

        int t1 = time_now();
        if (output_to_err_file == true){
            std::cerr << "initialize time: " << t1 - t0 << " seconds" << std::endl;
        }
    };

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // calculate the center and radius of the cell 
    void calculate_rc(std::vector<int>& atoms){

        Vector3d rc = Vector3d::Zero();
        Vector3d r = Vector3d::Zero();
        int size_atoms = atoms.size();
        
        MatrixXd crd_atom = MatrixXd::Zero(size_atoms,3);

        // add all atom coordinates to crd_atom
        for (int i = 0; i < size_atoms; i++){

            for (int j = 0; j < 3; j++){
                crd_atom(i,j) = t_crd(atoms[i]-1, j);
            }
        }

        // calculate center and radius of the cell
        for (int i = 0; i < 3; i++){

            double min_r = crd_atom.col(i).minCoeff();
            r(i) = abs(crd_atom.col(i).maxCoeff() - min_r) * 0.5;   // radius is half of abs(max - min)
            rc(i) = min_r + r(i);                                   // median of max and min
        }
        double max_r = r.maxCoeff();    // maximum radius among x, y, z

        // assign to structure Root_r
        root_r.rc = rc;
        root_r.r  = max_r;
    };

    // put new child cell to cells
    void add_child(int octant, int idx_current, std::vector<Cell>& cells){

        // insert new child cell to the end of cells
        cells.push_back(Cell());

        int idx_child = cells.size() - 1;

        // calculate center and radius of the new child cell
        cells[idx_child].r     = cells[idx_current].r * 0.5;
        cells[idx_child].rc(0) = cells[idx_current].rc(0) + cells[idx_child].r * ((octant & 1) * 2 - 1);
        cells[idx_child].rc(1) = cells[idx_current].rc(1) + cells[idx_child].r * ((octant & 2) - 1);
        cells[idx_child].rc(2) = cells[idx_current].rc(2) + cells[idx_child].r * ((octant & 4) / 2 - 1);

        cells[idx_child].parent = idx_current;

        // link new child cell to its parent cell
        cells[idx_current].child[octant] = idx_child;
        cells[idx_current].nchild        = (cells[idx_current].nchild | (1 << octant));
    };

    // split cell which has more than n_crit atoms (replace atoms in a cell to its child cells)
    void split_cell(int idx_current, std::vector<Cell>& cells) {

        std::queue<int> cell_queue;
        cell_queue.push(idx_current);
    
        // split cells until all child cells have n_crit or less atoms
        while (!(cell_queue.empty())) {
            int idx_curr_cell = cell_queue.front();
            cell_queue.pop();
    
            // reassign atoms in the cell to its child cells
            for (int i = 0; i < cells[idx_curr_cell].nleaf; i++) {

                int atom_current     = cells[idx_curr_cell].leaf[i];
                int idx_atom_current = atom_current - 1;

                int octant = (t_crd(idx_atom_current, 0) > cells[idx_curr_cell].rc(0)) + \
                             ((t_crd(idx_atom_current, 1) > cells[idx_curr_cell].rc(1)) << 1) + \
                             ((t_crd(idx_atom_current, 2) > cells[idx_curr_cell].rc(2)) << 2);   // determine octant
    
                // if the child cell does not exist, create it
                if (!(cells[idx_curr_cell].nchild & (1 << octant))) {
                    add_child(octant, idx_curr_cell, cells);
                }
                
                // put atom to the child cell
                int idx_child = cells[idx_curr_cell].child[octant];
                cells[idx_child].leaf.push_back(atom_current);
                cells[idx_child].nleaf += 1;

                if (cells[idx_child].nleaf > n_crit) {
                    cell_queue.push(idx_child);
                }
            }
        }
    };
    
    // create all_cells for all groups (main function to setup cells)
    std::vector<All_cells> setup_all_cells(){
        int t0 = time_now();

        all_cells = std::vector<All_cells>();
        int size_group = gname_iatoms_pairs.size();
        all_cells.reserve(size_group);

        for (int i = 0; i < size_group; i++){

            std::string group = gname_iatoms_pairs[i].first;
            std::vector<int> iatoms = gname_iatoms_pairs[i].second;
            int num_iatoms = iatoms.size();

            // set root cell
            All_cells all_cell;
            all_cell.group = group;
            all_cell.cells = std::vector<Cell>();

            std::vector<Cell>& cells = all_cell.cells;
            cells.push_back(Cell());
            calculate_rc(iatoms);
            cells[0].rc = root_r.rc;
            cells[0].r  = root_r.r;

            for (int j = 0; j < num_iatoms; j++){

                // start from the root cell
                int idx_current = 0;            // index of current cell
                int idx_atom = iatoms[j] - 1;

                while (cells[idx_current].nleaf > n_crit) {

                    // put atom to the current cell, then go to the child cell
                    cells[idx_current].leaf.push_back(iatoms[j]);
                    Vector3d crd_atom = t_crd.row(idx_atom);
                    cells[idx_current].nleaf += 1;

                    int octant = (crd_atom(0) > cells[idx_current].rc(0)) + \
                                 ((crd_atom(1) > cells[idx_current].rc(1)) << 1) + \
                                 ((crd_atom(2) > cells[idx_current].rc(2)) << 2);  // determine octant
             
                    // if the child cell does not exist, create it
                    if (!(cells[idx_current].nchild & (1 << octant))){
                        add_child(octant, idx_current, cells);
                    }
                    // go to the child cell
                    idx_current = cells[idx_current].child[octant];
                }
                // put atom to the current cell
                cells[idx_current].leaf.push_back(iatoms[j]);
                cells[idx_current].nleaf += 1;

                // if the current cell has more than n_crit atoms, split it
                if (cells[idx_current].nleaf > n_crit){
                    split_cell(idx_current, cells);
                }
            }
            all_cells.push_back(all_cell);
        }
        int t1 = time_now();
        if (output_to_err_file == true){
            std::cerr << "setup_all_cells time: " << t1 - t0 << " seconds" << std::endl;
        }
        return all_cells;
    };

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // set all_cells
    void get_all_cells(const std::vector<All_cells>& all_cells_input){

        all_cells = all_cells_input;  
    };    

    // calculate multipole/multipole_j of a cell
    void cal_multipole(VectorXd& multipole, MatrixXd& multipole_j, Vector3d& rc, std::vector<int>& atoms){
        
        int size_atoms = atoms.size();
        
        for (int i = 0; i < size_atoms; i++){

            int idx_atom = atoms[i] - 1;

            Vector3d crd_atom = t_crd.row(idx_atom);
            double dx = rc(0) - crd_atom(0);
            double dy = rc(1) - crd_atom(1);
            double dz = rc(2) - crd_atom(2);
            double qj = charges[idx_atom];

            double qjdx = qj * dx;
            double qjdy = qj * dy;
            double qjdz = qj * dz;
            
            VectorXd multipole_atom = VectorXd::Zero(10);

            // calculate multipole
            multipole_atom(0) = qj * 1.0;
            multipole_atom(1) = qjdx;
            multipole_atom(2) = qjdy;
            multipole_atom(3) = qjdz;
            multipole_atom(4) = qjdx * dx * 0.5;
            multipole_atom(5) = qjdy * dy * 0.5;
            multipole_atom(6) = qjdz * dz * 0.5;
            multipole_atom(7) = qjdx * dy;
            multipole_atom(8) = qjdy * dz;
            multipole_atom(9) = qjdz * dx;
            
            multipole += multipole_atom;

            // calculate multipole_j
            for (int j = 0; j < 3; j++){
                multipole_j.row(j) += multipole_atom.transpose() * crd_atom(j);
            }
        }
    };

    // calculate multipoles for all cells (only which do not have child cells)
    void cal_p(){
        int size_all_cells = all_cells.size();

        for (int i = 0; i < size_all_cells; i++){

            std::vector<Cell>& cells = all_cells[i].cells;
            int size_cells = cells.size();

            for (int j = 0; j < size_cells; j++){
               
                Cell& cell = cells[j];

                // skip if the cell has child cells
                if (cell.nchild != 0){
                    continue;
                }
                else {
                    cal_multipole(cell.multipole, cell.multipole_j, cell.rc, cell.leaf);
                }
            }
        }
    };

    // calculate multipoles of parent cells from child cells
    void cal_M2M(){

        int size_all_cells = all_cells.size();

        for (int i = 0; i < size_all_cells; i++){

            std::vector<Cell>& cells = all_cells[i].cells;
            int cells_size = cells.size();

            for (int i = 0; i < cells_size - 1; i++){

                int i_inv = cells_size - 1 - i;

                int idx_p = cells[i_inv].parent;    // index of parent cell
                int idx_c = i_inv;                  // index of child cell
            
                VectorXd& p_multipole = cells[idx_p].multipole;
                VectorXd& c_multipole = cells[idx_c].multipole;
                Vector3d& p_rc = cells[idx_p].rc;
                Vector3d& c_rc = cells[idx_c].rc;
                
                double dx = p_rc(0) - c_rc(0);
                double dy = p_rc(1) - c_rc(1);
                double dz = p_rc(2) - c_rc(2);
                double Mx = c_multipole(0) * dx;
                double My = c_multipole(0) * dy;
                double Mz = c_multipole(0) * dz;

                // calculate multipole of parent cell from child cell
                p_multipole(0) += c_multipole(0);
                p_multipole(1) += c_multipole(1) + Mx;
                p_multipole(2) += c_multipole(2) + My;
                p_multipole(3) += c_multipole(3) + Mz;
                p_multipole(4) += c_multipole(4) + dx * c_multipole(1) + 0.5 * Mx * dx;
                p_multipole(5) += c_multipole(5) + dy * c_multipole(2) + 0.5 * My * dy;
                p_multipole(6) += c_multipole(6) + dz * c_multipole(3) + 0.5 * Mz * dz;
                p_multipole(7) += c_multipole(7) + dy * c_multipole(1) + dx * c_multipole(2) + Mx * dy;
                p_multipole(8) += c_multipole(8) + dz * c_multipole(2) + dy * c_multipole(3) + My * dz;
                p_multipole(9) += c_multipole(9) + dx * c_multipole(3) + dz * c_multipole(1) + Mz * dx;

                MatrixXd& p_multipole_j = cells[idx_p].multipole_j;
                MatrixXd& c_multipole_j = cells[idx_c].multipole_j;

                // calculate multipole_j of parent cell from child cell
                for (int j = 0; j < 3; j++){

                    double Mjx = c_multipole_j(j, 0) * dx;
                    double Mjy = c_multipole_j(j, 0) * dy;
                    double Mjz = c_multipole_j(j, 0) * dz;

                    p_multipole_j(j, 0) += c_multipole_j(j, 0);
                    p_multipole_j(j, 1) += c_multipole_j(j, 1) + Mjx;
                    p_multipole_j(j, 2) += c_multipole_j(j, 2) + Mjy;
                    p_multipole_j(j, 3) += c_multipole_j(j, 3) + Mjz;
                    p_multipole_j(j, 4) += c_multipole_j(j, 4) + dx * c_multipole_j(j, 1) + 0.5 * Mjx * dx;
                    p_multipole_j(j, 5) += c_multipole_j(j, 5) + dy * c_multipole_j(j, 2) + 0.5 * Mjy * dy;
                    p_multipole_j(j, 6) += c_multipole_j(j, 6) + dz * c_multipole_j(j, 3) + 0.5 * Mjz * dz;
                    p_multipole_j(j, 7) += c_multipole_j(j, 7) + dy * c_multipole_j(j, 1) + dx * c_multipole_j(j, 2) + Mjx * dy;
                    p_multipole_j(j, 8) += c_multipole_j(j, 8) + dz * c_multipole_j(j, 2) + dy * c_multipole_j(j, 3) + Mjy * dz;
                    p_multipole_j(j, 9) += c_multipole_j(j, 9) + dx * c_multipole_j(j, 3) + dz * c_multipole_j(j, 1) + Mjz * dx;
                }
            }
        }
    };

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // determine near and far cells for each source atom
    void determine_near_far(int num_source, int idx_p, std::vector<Cell>& cells){

        int      idx_source = num_source - 1;
        Vector3d crd_source = t_crd.row(idx_source);

        std::queue<int> cell_queue;
        cell_queue.push(idx_p);

        while (!cell_queue.empty()){
            int idx_p = cell_queue.front();                     // index of current parent cell
            cell_queue.pop();
        
            if (cells[idx_p].nleaf > n_crit){

                for (int octant = 0; octant < 8; octant++){
                    
                    if (cells[idx_p].nchild & (1 << octant)) {
                        
                        int       idx_c = cells[idx_p].child[octant]; // index of current child cell
                        Vector3d& rc    = cells[idx_c].rc;
                        
                        double dx = crd_source(0) - rc(0);
                        double dy = crd_source(1) - rc(1);
                        double dz = crd_source(2) - rc(2);
                        double r  = sqrt(dx * dx + dy * dy + dz * dz);

                        double cutoff  = 9.0;
                        double min_len = r - 1.7320508 * cells[idx_c].r;        // the minimum length of particle(source) - particle(target in cell[c]) 

                        if (min_len < cutoff || cells[idx_c].r > theta * r){

                            cell_queue.push(idx_c);
                        }
                        else{
                            pairs_far[num_source].second.push_back(idx_c);
                        }      
                    }
                }
            }
            else {
                pairs_near[num_source].second.push_back(idx_p);
            }
        }
    };

    // calculate force and flux for near atom-atom pairs
    void cal_force_near(std::vector<Cell>& cells){
        int t0 = time_now();

        int size_source = pairs_near.size();
        for (int i = 0; i < size_source; i++){

            int num_source = pairs_near[i].first;
            std::vector<int>& pairs = pairs_near[i].second;

            int      idx_source = num_source - 1;
            Vector3d crd_source = t_crd.row(idx_source);
            Vector3d vi         = t_vel.row(idx_source);

            int size_pairs = pairs.size();
            for (int j = 0; j < size_pairs; j++){

                int idx_cell = pairs[j];
                std::vector<int>& atom_pairs = cells[idx_cell].leaf;
                int num_target = atom_pairs.size();
                
                for (int k = 0; k < num_target; k++){

                    int idx_target = atom_pairs[k] - 1;

                    if (idx_target == idx_source){
                        continue;
                    }

                    Vector3d crd_target = t_crd.row(idx_target);
                    Vector3d rij = crd_source - crd_target;
                
                    double r = sqrt(rij.dot(rij));

                    double inv_r   = 1.0 / r;
                    double coeff   = 332.05221729;
                    double inv_r3  = inv_r * inv_r * inv_r;             // 1/r^3
                    double inv_r8  = inv_r3 * inv_r3 * inv_r * inv_r;   // 1/r^8
                    double inv_r14 = inv_r8 * inv_r3 * inv_r3;          // 1/r^14

                    double charge_i = charges[idx_source];
                    double charge_j = charges[idx_target];
                    double coulomb  = coeff * charge_i * charge_j * inv_r3;

                    int type_source = atom_types[idx_source] - 1;
                    int type_target = atom_types[idx_target] - 1;

                    double c6  = c6s(type_source, type_target);
                    double c12 = c12s(type_source, type_target);
                    double vdw = 12.0 * c12 * inv_r14 - 6.0 * c6 * inv_r8;
                
                    Vector3d fij = rij * (coulomb + vdw); // electrostatic + vdw force

                    cal_flux_near(fij, vi, rij, idx_source, idx_target);
                }
            }
        }
        int t1 = time_now();
        if (output_to_err_file == true){
            std::cerr << "cal_force_near time: " << t1 - t0 << std::endl;
        }
    };

    // calculate force and flux for far atom-cell pairs
    void cal_force_far(std::vector<Cell>& cells){

        int t0 = time_now();

        int size_source = pairs_far.size();
        Matrix<double, 10, 1> bJx, bJy, bJz;
        Matrix<double, 10, 3> bJ;

        for (int i = 0; i < size_source; i++){

            int num_source = pairs_far[i].first;
            int idx_source = num_source - 1;

            Vector3d crd_source = t_crd.row(idx_source);
            Vector3d vi         = t_vel.row(idx_source);
            int size_pairs      = pairs_far[i].second.size();

            for (int j = 0; j < size_pairs; j++){

                int   idx_cell = pairs_far[i].second[j];
                Cell& cell     = cells[idx_cell];

                Vector3d rc = cell.rc;
                double   dx = crd_source(0) - rc(0);
                double   dy = crd_source(1) - rc(1);
                double   dz = crd_source(2) - rc(2);
                double   r  = sqrt(dx * dx + dy * dy + dz * dz);

                bJx.setZero();
                bJy.setZero();
                bJz.setZero();
                bJ.setZero();

                double inv_r = 1.0 / r;                     // 1/r
                double r2 = inv_r * inv_r;                  // 1/r^2
                double r3 = r2 * inv_r;                     // 1/r^3
                double r5 = r3 * r2;                        // 1/r^5
                double r7 = r5 * r2;                        // 1/r^7

                double dx2 = dx * dx;                       // dx^2
                double dy2 = dy * dy;                       // dy^2
                double dz2 = dz * dz;                       // dz^2

                double dxdy = dx * dy;                      // dxdy
                double dydz = dy * dz;                      // dydz
                double dzdx = dz * dx;                      // dzdx

                double dxr5 = 3 * dx * r5;                  // 3dx/r^5
                double dyr5 = 3 * dy * r5;                  // 3dy/r^5
                double dzr5 = 3 * dz * r5;                  // 3dx/r^5

                double dxdydz = 15 * dxdy * dz * r7;        // 15dxdydz/r^7
                double dx2dy  = 15 * dx2 * dy * r7;         // 15dx^2dy/r^7
                double dy2dz  = 15 * dy2 * dz * r7;         // 15dy^2dz/r^7
                double dz2dx  = 15 * dz2 * dx * r7;         // 15dz^2dx/r^7
                double dy2dx  = 15 * dy2 * dx * r7;         // 15dy^2dx/r^7
                double dz2dy  = 15 * dz2 * dy * r7;         // 15dz^2dy/r^7
                double dx2dz  = 15 * dx2 * dz * r7;         // 15dx^2dz/r^7

                // calculate bJx
                bJx(0) = -dx * r3;                          // -dx/r^3
                bJx(1) = -r3 + 3 * dx2 * r5;                // -1/r^3  + 3dx^2/r^5
                bJx(2) = 3 * dxdy * r5;                     // 0       + 3dydx/r^5
                bJx(3) = 3 * dzdx * r5;                     // 0       + 3dzdx/r^5
                bJx(4) = 3 * dxr5 - 15 * dx2 * dx * r7;     // 9dx/r^5 - 15dx^3/r^7
                bJx(5) =     dxr5 - dy2dx;                  // 3dx/r^5 - 15dy^2dx/r^7
                bJx(6) =     dxr5 - dz2dx;                  // 3dx/r^5 - 15dz^2dx/r^7
                bJx(7) =     dyr5 - dx2dy;                  // 3dy/r^5 - 15dx^2dy/r^7
                bJx(8) =          - dxdydz;                 // 0       - 15dydzdx/r^7
                bJx(9) =     dzr5 - dx2dz;                  // 3dz/r^5 - 15dx^2dz/r^7

                // calculate bJy
                bJy(0) = -dy * r3;                          // -dy/r^3
                bJy(1) = bJx(2);                            // 0       + 3dxdy/r^5
                bJy(2) = -r3 + 3 * dy2 * r5;                // -1/r^3  + 3dy^2/r^5
                bJy(3) = 3 * dydz * r5;                     // 0       + 3dzdy/r^5
                bJy(4) =     dyr5 - dx2dy;                  // 3dy/r^5 - 15dx^2dy/r^7
                bJy(5) = 3 * dyr5 - 15 * dy2 * dy * r7;     // 9dy/r^5 - 15dy^3/r^7
                bJy(6) =     dyr5 - dz2dy;                  // 3dy/r^5 - 15dz^2dy/r^7
                bJy(7) =     dxr5 - dy2dx;                  // 3dx/r^5 - 15dxdy^2/r^7
                bJy(8) =     dzr5 - dy2dz;                  // 3dz/r^5 - 15dydzdy/r^7
                bJy(9) = bJx(8);                            // 0       - 15dzdxdy/r^7

                // calculate bJz
                bJz(0) = -dz * r3;                          // -dz/r^3
                bJz(1) = bJx(3);                            // 0       + 3dxdz/r^5
                bJz(2) = bJy(3);                            // 0       + 3dydz/r^5
                bJz(3) = -r3 + 3 * dz2 * r5;                // -1/r^3  + 3dz^2/r^5
                bJz(4) =     dzr5 - dx2dz;                  // 3dz/r^5 - 15dx^2dz/r^7
                bJz(5) =     dzr5 - dy2dz;                  // 3dz/r^5 - 15dy^2dz/r^7
                bJz(6) = 3 * dzr5 - 15 * dz2 * dz * r7;     // 9dz/r^5 - 15dz^3/r^7
                bJz(7) = bJx(8);                            // 0       - 15dxdydz/r^7
                bJz(8) =     dyr5 - dz2dy;                  // 3dy/r^5 - 15dydz^2/r^7
                bJz(9) =     dxr5 - dz2dx;                  // 3dx/r^5 - 15dzdxdz/r^7

                // calculate multipole
                VectorXd& multipole = cell.multipole;
                double charge = charges[idx_source];
                double coeff  = 332.05221729;
                VectorXd multi = coeff * multipole * charge;

                double   fx = multi.dot(bJx);
                double   fy = multi.dot(bJy);
                double   fz = multi.dot(bJz);
                Vector3d f  = Vector3d(-fx, -fy, -fz);

                MatrixXd& multipole_j = cell.multipole_j;
                bJ.col(0) = bJx;
                bJ.col(1) = bJy;
                bJ.col(2) = bJz;
                Matrix3d multi_j = -multipole_j * bJ * coeff * charge;

                int igrp = iatom_to_igroup(idx_source);
                int Jgrp = iatom_to_igroup(cell.leaf[0] - 1);
                Vector3d vi = t_vel.row(idx_source);

                cal_flux_far(f, vi, crd_source, multi_j, igrp, Jgrp);
            }
        }
        int t1 = time_now();
        if (output_to_err_file == true){
            std::cerr << "cal_force_far time: " << t1 - t0 << std::endl;
        }
    };

    // calculate heat flux for near atom-atom pairs
    void cal_hflux_near(Vector3d fij, Vector3d vi, Vector3d rij, int idx_source, int idx_target){
        
        Vector3d h_ij = rij * (fij.dot(vi)) * 0.5;

        int igrp = iatom_to_igroup(idx_source);
        int jgrp = iatom_to_igroup(idx_target);

        if (igrp != jgrp){
            hflux_ij_near[igrp-1][jgrp-1] += h_ij;
            hflux_ij_near[jgrp-1][igrp-1] += h_ij;
        }
        else {
            hflux_ij_near[igrp-1][jgrp-1] += h_ij;
        }
        count_near += 1;
    };

    // calculate heat flux for far atom-cell pairs
    void cal_hflux_far(Vector3d fiJ, Vector3d vi, Vector3d r, Matrix3d multi_j, int igrp, int Jgrp){
        
        Vector3d h_ij = r * (fiJ.dot(vi)) * 0.5 - multi_j * vi * 0.5;

        if (igrp != Jgrp){
            hflux_ij_far[igrp-1][Jgrp-1] += h_ij;
            hflux_ij_far[Jgrp-1][igrp-1] += h_ij;
        }
        else{
            hflux_ij_far[igrp-1][Jgrp-1] += h_ij;
        }
        count_far += 1;
    };

    // calculate energy flux for near atom-atom pairs
    void cal_eflux_near(Vector3d fij, Vector3d vi, Vector3d rij, int idx_source, int idx_target){

        double e_ij = fij.dot(vi) * 0.5;

        int igrp = iatom_to_igroup(idx_source);
        int jgrp = iatom_to_igroup(idx_target);

        if (igrp != jgrp){
            eflux_ij_near(igrp-1, jgrp-1) += e_ij;
            eflux_ij_near(jgrp-1, igrp-1) -= e_ij;
        }
        else if (idx_source < idx_target) {
            eflux_ij_near(igrp-1, jgrp-1) += e_ij;
        }
        else {
            eflux_ij_near(igrp-1, jgrp-1) -= e_ij;
        }
        count_near += 1;
    };

    // calculate energy flux for far atom-cell pairs
    void cal_eflux_far(Vector3d fiJ, Vector3d vi, Vector3d r, Matrix3d multi_j, int igrp, int Jgrp){
        
        double e_ij = fiJ.dot(vi) * 0.5;

        if (igrp != Jgrp){
            eflux_ij_far(igrp-1, Jgrp-1) += e_ij;
            eflux_ij_far(Jgrp-1, igrp-1) -= e_ij;
        }
        else{
            eflux_ij_far(igrp-1, Jgrp-1) += e_ij;
        }
        count_far += 1;
    };

    // get atoms vector of a group
    std::vector<int> get_atoms(const std::string& source, const std::vector<std::pair<std::string, std::vector<int>>>& pairs) {

        std::vector<int> atoms;
        for (const auto& pair : pairs) {

            if (pair.first == source) {
                atoms = pair.second;
                break;
            }
        }
        return atoms;
    };

    // get cells vector of a group
    std::vector<Cell> get_cells(const std::string& source, const std::vector<All_cells>& all_cells) {
        
        std::vector<Cell> cells;
        for (const auto& all_cell : all_cells) {

            if (all_cell.group == source) {
                cells =  all_cell.cells;
                break;
            }
        }
        return cells;
    };

    // calculate forces and fluxes for all required group pairs
    void cal_force_and_flux(){

        int size_table = gpair_table.size();
        for (int i = 0; i < size_table; i++){

            std::string source = gpair_table[i].first;
            std::vector<std::string> targets = gpair_table[i].second;
            int size_targets = targets.size();
        
            std::vector<int> source_atoms = get_atoms(source, gname_iatoms_pairs);
            int source_size = source_atoms.size();

            for (int j = 0; j < size_targets; j++){

                pairs_near.clear();
                pairs_far.clear();
                pairs_near.resize(source_size);
                pairs_far.resize(source_size);
                std::string target = targets[j];
                std::vector<Cell> cells = get_cells(target, all_cells);

                for (int k = 0; k < source_size; k++){

                    int num_source = source_atoms[k];
                    pairs_near[k].first = num_source;
                    pairs_far[k].first = num_source;
                    determine_near_far(num_source, 0, cells);
                }
                cal_force_near(cells);
                cal_force_far(cells);
            }

            std::vector<Cell> cells_source = get_cells(source, all_cells);
            for (int l = 0; l < size_targets; l++){

                std::string target = targets[l];
                if (target == source){
                    continue;
                }
                std::vector<int> target_atoms = get_atoms(target, gname_iatoms_pairs);
                int target_size = target_atoms.size();
                pairs_near.clear();
                pairs_far.clear();
                pairs_near.resize(target_size);
                pairs_far.resize(target_size);

                for (int m = 0; m < target_size; m++){
                    int num_target = target_atoms[m];
                    determine_near_far(num_target, 0, cells_source);
                }
                cal_force_near(cells_source);
                cal_force_far(cells_source);
            }
        }
    };

    // calculate forces and fluxes for bonded atom pairs (to exclude from calculation)
    void cal_bonded_flux(){

        int size_bonded = bonded_pairs.size();

        for (int i = 0; i < size_bonded; i++){

            int num_source = bonded_pairs[i].first;
            int num_target = bonded_pairs[i].second;
            int idx_source = num_source - 1;
            int idx_target = num_target - 1;

            Vector3d crd_source = t_crd.row(idx_source);
            Vector3d crd_target = t_crd.row(idx_target);

            Vector3d rij = crd_source - crd_target;
            Vector3d vi  = t_vel.row(idx_source);
            Vector3d vj  = t_vel.row(idx_target);
            Vector3d vij = vi + vj;

            double r       = sqrt(rij.dot(rij));
            double inv_r   = 1.0 / r;
            double inv_r3  = inv_r * inv_r * inv_r;              // 1/r^3
            double inv_r8  = inv_r3 * inv_r3 * inv_r * inv_r;    // 1/r^8
            double inv_r14  = inv_r8 * inv_r3 * inv_r3;          // 1/r^14
            
            double coeff    = 332.05221729;
            double charge_i = charges[idx_source];
            double charge_j = charges[idx_target];
            double coulomb  = coeff * charge_i * charge_j * inv_r3;

            int type_source = atom_types[idx_source] - 1;
            int type_target = atom_types[idx_target] - 1;


            double c6  = c6s(type_source, type_target);
            double c12 = c12s(type_source, type_target);
            double vdw = 12.0 * c12 * inv_r14 - 6.0 * c6 * inv_r8;

            Vector3d fij = rij * (coulomb + vdw);

            int igrp = iatom_to_igroup(idx_source);
            int jgrp = iatom_to_igroup(idx_target);

            cal_flux_near_bonded(fij, vij, rij, igrp, jgrp);
        }
        if (output_to_err_file == true){
            std::cerr << "count_atom_bonded: " << size_bonded << std::endl;
        }
    };

    // calculate heat flux for bonded atom pairs
    void cal_hflux_near_bonded(Vector3d fij, Vector3d vij, Vector3d rij, int igrp, int jgrp){
        
        Vector3d h_ij = rij * (fij.dot(vij)) * 0.5;

        if (igrp != jgrp){
            hflux_ij_near[igrp-1][jgrp-1] -= h_ij;
            hflux_ij_near[jgrp-1][igrp-1] -= h_ij;
        }
        else{
            hflux_ij_near[igrp-1][jgrp-1] -= h_ij;
        }
    };

    // calculate energy flux for bonded atom pairs
    void cal_eflux_near_bonded(Vector3d fij, Vector3d vij, Vector3d rij, int igrp, int jgrp){

        double e_ij = fij.dot(vij) * 0.5;
        
        if (igrp != jgrp){
            eflux_ij_near(igrp-1, jgrp-1) -= e_ij;
            eflux_ij_near(jgrp-1, igrp-1) += e_ij;
        }
        else{
            eflux_ij_near(igrp-1, jgrp-1) -= e_ij;
        }
    };

    // add far and near flux
    void add_flux(){

        if (flag_heat == true){
            for (int i = 0; i < ngrp; i++){
                for (int j = 0; j < ngrp; j++){
                    hflux_ij[i][j] = hflux_ij_far[i][j] + hflux_ij_near[i][j];
                }
            }
        }
        else if (flag_energy == true){
            eflux_ij = eflux_ij_far + eflux_ij_near;
        }
    };

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // main function to calculate coulomb flux using fmm
    void cal_coulomb_flux_fmm(const std::vector<All_cells>& all_cells) {
        int t0 = time_now();
        get_all_cells(all_cells);
        int t1 = time_now();
        cal_p();
        int t2 = time_now();
        cal_M2M();
        int t3 = time_now();
        cal_force_and_flux();
        int t4 = time_now();
        cal_bonded_flux();
        int t5 = time_now();
        add_flux();
        int t6 = time_now();
        if (output_to_err_file == true){        
            std::cerr << "get_all_cells time: " << t1 - t0 << std::endl;
            std::cerr << "cal_p time: " << t2 - t1 << std::endl;
            std::cerr << "cal_M2M time: " << t3 - t2 << std::endl;
            std::cerr << "cal_force_and_flux time: " << t4 - t3 << std::endl;
            std::cerr << "cal_bonded_flux time: " << t5 - t4 << std::endl;
            std::cerr << "add_flux time: " << t6 - t5 << std::endl; 
            std::cerr << "total time(cal_p -> add_flux): " << t6 - t1 << std::endl;
            std::cerr << "count_far: " << count_far << std::endl;
            std::cerr << "count_near: " << count_near << std::endl;
            std::cerr << "ncrit: " << n_crit << std::endl;
            std::cerr << "theta: " << theta << std::endl;
            std::cerr << "   " << std::endl;
        }
    };

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // get current time in microseconds
    int time_now(){
 
        auto start = std::chrono::system_clock::now();
        auto now = std::chrono::duration_cast<std::chrono::microseconds>(start.time_since_epoch()).count();
        return now;
    };
};

// Convert hflux_ij to numpy array
py::array_t<double> hflux_to_numpy(const std::vector<std::vector<Vector3d>>& hflux_ij) {
    size_t ngrp_i = hflux_ij.size();
    size_t ngrp_j = hflux_ij.empty() ? 0 : hflux_ij[0].size();

    // shape: (ngrp_i, ngrp_j, 3)
    std::vector<ssize_t> shape = { static_cast<ssize_t>(ngrp_i),
                                   static_cast<ssize_t>(ngrp_j),
                                   3 };
    std::vector<ssize_t> strides = {
        static_cast<ssize_t>(ngrp_j * 3 * sizeof(double)),
        static_cast<ssize_t>(3 * sizeof(double)),
        static_cast<ssize_t>(sizeof(double))
    };

    // Prepare the data (flat array)
    std::vector<double> buffer;
    buffer.reserve(ngrp_i * ngrp_j * 3);

    for (size_t i = 0; i < ngrp_i; ++i) {
        for (size_t j = 0; j < ngrp_j; ++j) {
            const Vector3d& v = hflux_ij[i][j];
            buffer.push_back(v(0));
            buffer.push_back(v(1));
            buffer.push_back(v(2));
        }
    }

    // Return as numpy array (make a copy to ensure safety)
    return py::array(py::buffer_info(
        buffer.data(),
        sizeof(double),
        py::format_descriptor<double>::format(),
        3,
        shape,
        strides
    )).attr("copy")();
};

// debug function(print all cells)
void print_all_cells(){
    int size_all_cells = all_cells.size();
    for (int i = 0; i < size_all_cells; i++){
        std::cerr << "group: " << all_cells[i].group << std::endl;
        std::vector<Cell>& cells = all_cells[i].cells;
        int size_cells = cells.size();
        for (int j = 0; j < size_cells; j++){
            std::cerr << "cell: " << j << " nleaf: " << cells[j].nleaf << " leaves: [" ;

            for (size_t k = 0; k < cells[j].leaf.size(); ++k) {
                std::cerr << cells[j].leaf[k];
                if (k < cells[j].leaf.size() - 1) {
                    std::cerr << ", ";
                }
            };
            std::cerr << "]" << std::endl;

            std::cerr << "nchild: " << cells[j].nchild << " child: [" ;
            for (size_t k = 0; k < cells[j].child.size(); ++k) {
                std::cerr << cells[j].child[k];
                if (k < cells[j].child.size() - 1) {
                    std::cerr << ", ";
                }
            }
            std::cerr << "]" << std::endl;
            std::cerr << "parent: " << cells[j].parent << std::endl;
            std::cerr << "rc: " << cells[j].rc.transpose() << std::endl;
            std::cerr << "r: " << cells[j].r << std::endl;
            if (cells[j].multipole(0) != 0){
                std::cerr << "multipole: " << cells[j].multipole.transpose() << std::endl;
            }
            if (cells[j].nchild != 0 && cells[j].multipole[0] != 0){
                std::cerr << "error: multipole exists: [";

                for (size_t k = 0; k < cells[j].leaf.size(); ++k) {
                    std::cerr << cells[j].leaf[k];
                    if (k < cells[j].leaf.size() - 1) {
                        std::cerr << ", ";
                    }
                };
                std::cerr << "]" << std::endl;
            }
            std::cerr << "   " << std::endl;
        }
    }
};


PYBIND11_MODULE(lib_flux_fmm, m){
    
    py::class_<cal_fmm>(m, "cal_fmm")
        .def(py::init<>())
        .def("setup", &cal_fmm::setup)
        .def("set_flux", &cal_fmm::set_flux)
        .def("initialize", &cal_fmm::initialize)
        .def("setup_all_cells", &cal_fmm::setup_all_cells)
        .def("cal_coulomb_flux_fmm", &cal_fmm::cal_coulomb_flux_fmm)
        .def("check_cells", &cal_fmm::check_cells)
        .def_readwrite("natom", &cal_fmm::natom)
        .def_readwrite("n_crit", &cal_fmm::n_crit)
        .def_readwrite("theta", &cal_fmm::theta)
        .def_readwrite("charges", &cal_fmm::charges)
        .def_readwrite("t_crd", &cal_fmm::t_crd)
        .def_readwrite("t_vel", &cal_fmm::t_vel)
        .def_readwrite("bonded_pairs", &cal_fmm::bonded_pairs)
        .def_readwrite("gpair_table", &cal_fmm::gpair_table)
        .def_readwrite("gname_iatoms_pairs", &cal_fmm::gname_iatoms_pairs)
        .def_readwrite("iatom_to_igroup", &cal_fmm::iatom_to_igroup)
        .def_readwrite("all_cells", &cal_fmm::all_cells)
        .def_readwrite("hflux_ij", &cal_fmm::hflux_ij)
        .def_readwrite("eflux_ij", &cal_fmm::eflux_ij)
        ;
    
    
    py::class_<cal_fmm::Cell>(m, "Cell")
        .def(py::init<>())
        .def_readwrite("nleaf", &cal_fmm::Cell::nleaf)
        .def_readwrite("leaf", &cal_fmm::Cell::leaf)
        .def_readwrite("nchild", &cal_fmm::Cell::nchild)
        .def_readwrite("child", &cal_fmm::Cell::child)
        .def_readwrite("parent", &cal_fmm::Cell::parent)
        .def_readwrite("rc", &cal_fmm::Cell::rc)
        .def_readwrite("r", &cal_fmm::Cell::r)
        .def_readwrite("multipole", &cal_fmm::Cell::multipole)
        ;
    
    py::class_<cal_fmm::All_cells>(m, "All_cells")
        .def(py::init<>())
        .def_readwrite("group", &cal_fmm::All_cells::group)
        .def_readwrite("cells", &cal_fmm::All_cells::cells)
        ;

    m.def("hflux_to_numpy", &hflux_to_numpy, "Convert hflux_ij to numpy array");
    m.def("print_all_cells", &print_all_cells, "Print all cells for debug");
}