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
    
public:
    // public variables
    int natom;
    int ngrp;
    int n_crit;
    float theta;
    std::vector<double> charges;
    MatrixXd t_crd;
    MatrixXd t_vel;

    std::vector<std::vector<std::tuple<int, int, int>>> interact_table;
    std::vector<std::pair<std::string, std::vector<std::string>>> gpair_table;
    std::vector<std::pair<std::string, std::vector<int>>> gname_iatoms_pairs;
    VectorXi iatom_to_igroup;

    std::vector<std::vector<Vector3d>> hflux_ij;
    MatrixXd eflux_ij;
    bool flag_heat;
    bool flag_energy;

    int count_atom;
    int count_cell;

    std::function<void(Vector3d, Vector3d, Vector3d, int, int)> cal_flux_atomwise;
    std::function<void(Vector3d, Vector3d, Vector3d, Matrix3d, int, int)> cal_flux_cellwise;

    // constructor
    cal_fmm():
        natom(0),
        ngrp(0),
        n_crit(0),
        theta(0.0),
        charges(std::vector<double>()),
        t_crd(MatrixXd::Zero(natom, 3)),
        t_vel(MatrixXd::Zero(natom, 3)),
        interact_table(std::vector<std::vector<std::tuple<int, int, int>>>()),
        gpair_table(std::vector<std::pair<std::string, std::vector<std::string>>>()),
        gname_iatoms_pairs(std::vector<std::pair<std::string, std::vector<int>>>()),
        iatom_to_igroup(VectorXi::Zero(0)),
        hflux_ij(std::vector<std::vector<Vector3d>>()),
        eflux_ij(MatrixXd::Zero(0, 0)),
        flag_heat(false),
        flag_energy(false),
        count_atom(0),
        count_cell(0)
    {};

    void setup(const int& input_natom, const int& input_n_crit, const float& input_theta, const std::vector<double>& input_charges, \
        const std::vector<std::vector<std::tuple<int, int, int>>>& input_interact_table, \
        const std::vector<std::pair<std::string, std::vector<int>>>& input_gname_iatoms_pairs, \
        const std::vector<std::pair<std::string, std::vector<std::string>>>& input_gpair_table, \
        const VectorXi& input_iatom_to_igroup){
        int t0 = time_now();
        natom  = input_natom;
        n_crit = input_n_crit;
        theta  = input_theta;
        charges = input_charges;
        interact_table = input_interact_table;
        gname_iatoms_pairs = input_gname_iatoms_pairs;
        gpair_table = input_gpair_table;
        iatom_to_igroup = input_iatom_to_igroup;

        ngrp = iatom_to_igroup.maxCoeff();
        int t1 = time_now();
        std::cerr << "setup time: " << t1 - t0 << " seconds" << std::endl;
    };

    // get flux type
    void set_flux(const std::string& flux_type){

        if (flux_type == "heat"){
            flag_heat = true;
            hflux_ij.resize(ngrp, std::vector<Vector3d>(ngrp, Vector3d::Zero()));
            cal_flux_cellwise = [=](const Vector3d fiJ, const Vector3d vi, const Vector3d r, const Matrix3d pot_j, const int igrp, const int Jgrp) {
                cal_hflux_cellwise(fiJ, vi, r, pot_j, igrp, Jgrp);
            };
            cal_flux_atomwise = [=](const Vector3d fij, const Vector3d vi, const Vector3d rij, const int igrp, const int jgrp) {
                cal_hflux_atomwise(fij, vi, rij, igrp, jgrp);
            };
        }
        else if (flux_type == "energy"){
            flag_energy = true;
            eflux_ij = MatrixXd::Zero(ngrp, ngrp);
            cal_flux_cellwise = [=](const Vector3d fiJ, const Vector3d vi, const Vector3d r, const Matrix3d pot_j, const int igrp, const int Jgrp) {
                cal_eflux_cellwise(fiJ, vi, r, pot_j, igrp, Jgrp);
            };
            cal_flux_atomwise = [=](const Vector3d fij, const Vector3d vi, const Vector3d rij, const int igrp, const int jgrp) {
                cal_eflux_atomwise(fij, vi, rij, igrp, jgrp);
            };
        }
        else{
            std::cerr << "Invalid flux type" << std::endl;
            std::exit(1);
        }
    };

    // read trajectory
    void initialize(const MatrixXd& crd, const MatrixXd& vel){
        int t0 = time_now();
        t_crd = MatrixXd::Zero(natom, 3);
        t_vel = MatrixXd::Zero(natom, 3);
        t_crd = crd;
        t_vel = vel;
        count_atom = 0;
        count_cell = 0;
        if (flag_heat == true){
            hflux_ij.clear();
            hflux_ij.resize(ngrp, std::vector<Vector3d>(ngrp, Vector3d::Zero()));
        }
        else if (flag_energy == true){
            eflux_ij = MatrixXd::Zero(ngrp, ngrp);
        }
        int t1 = time_now();
        std::cerr << "initialize time: " << t1 - t0 << " seconds" << std::endl;

    };

////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Cell {
        int               nleaf;
        std::vector<int>  leaf;
        int               nchild;
        std::vector<int>  child;
        int               parent;
        Vector3d          rc;
        double            r;
        Vector3d          rmc;
        VectorXd          multipole;
        MatrixXd          multipole_j;

        Cell(int n_crit):
            nleaf(0),                               // number of atoms(leaf) in the cell
            leaf(std::vector<int>(0)),              // index of atoms in the cell
            nchild(0),                              // number of child cells
            child(std::vector<int>(8, 0)),          // index of 8 child cells
            parent(0),                              // index of parent cell
            rc(Vector3d::Zero()),                   // center of the cell
            r(0.0),                                 // radius of the cell
            rmc(Vector3d::Zero()),                  // center of mass of the cell
            multipole(VectorXd::Zero(10)),          // 10 multipoles
            multipole_j(MatrixXd::Zero(3, 10))
        {}
    };

    struct All_cells {
        std::string group;
        std::vector<Cell> cells;

        All_cells():
            group(""),
            cells()
        {}
    };
    std::vector<All_cells> all_cells;

////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct rc_max_r{
        Vector3d rc;
        double r;

        rc_max_r():
            rc(Vector3d::Zero()),
            r(0.0)
        {}
    };
    rc_max_r result_root;

    // calculate the center and radius of the cell 
    void calculate_rc(std::vector<int>& atoms){

        Vector3d rc = Vector3d::Zero();
        Vector3d r = Vector3d::Zero();
        int size_atoms = atoms.size();
        
        MatrixXd crd_atom = MatrixXd::Zero(size_atoms,3);

        for (int i = 0; i < size_atoms; i++){

            for (int j = 0; j < 3; j++){
                crd_atom(i,j) = t_crd(atoms[i]-1, j);
            }
        }

        for (int i = 0; i < 3; i++){

            double min_r = crd_atom.col(i).minCoeff();
            r(i) = abs(crd_atom.col(i).maxCoeff() - min_r) * 0.5;
            rc(i) = min_r + r(i);
        }
        double max_r = r.maxCoeff();

        result_root.rc = rc;
        result_root.r = max_r;
    };

    void add_child(int octant, int current_cell, std::vector<Cell>& cells){

        cells.push_back(Cell(n_crit));

        int new_child = cells.size() - 1;

        // std::cerr << "new child: " << new_child << std::endl;

        cells[new_child].r = cells[current_cell].r * 0.5;
        cells[new_child].rc(0) = cells[current_cell].rc(0) + cells[new_child].r * ((octant & 1) * 2 - 1);
        cells[new_child].rc(1) = cells[current_cell].rc(1) + cells[new_child].r * ((octant & 2) - 1);
        cells[new_child].rc(2) = cells[current_cell].rc(2) + cells[new_child].r * ((octant & 4) / 2 - 1);

        cells[new_child].parent = current_cell;
        cells[current_cell].child[octant] = new_child;
        cells[current_cell].nchild = (cells[current_cell].nchild | (1 << octant));
    };

    void split_cell(int current_cell, std::vector<Cell>& cells) {

        std::queue<int> cell_queue;
        cell_queue.push(current_cell);
    
        while (!(cell_queue.empty())) {
            int cell_index = cell_queue.front();
            cell_queue.pop();
    
            for (int i = 0; i < cells[cell_index].nleaf; i++) {
                int atom_current = cells[cell_index].leaf[i];
                int index_atom_current = atom_current - 1;
                int octant = (t_crd(index_atom_current, 0) > cells[cell_index].rc(0)) + \
                             ((t_crd(index_atom_current, 1) > cells[cell_index].rc(1)) << 1) + \
                             ((t_crd(index_atom_current, 2) > cells[cell_index].rc(2)) << 2);
    
                if (!(cells[cell_index].nchild & (1 << octant))) {
                    add_child(octant, cell_index, cells);
                }
    
                int child_cell = cells[cell_index].child[octant];
                cells[child_cell].leaf.push_back(atom_current);
                cells[child_cell].nleaf += 1;
                cells[child_cell].rmc += t_crd.row(index_atom_current);

                if (cells[child_cell].nleaf > n_crit) {
                    cell_queue.push(child_cell);
                }
            }
        }
    };
        
    std::vector<All_cells> setup_all_cells(){
        int t0 = time_now();

        all_cells = std::vector<All_cells>();
        int group_size = gname_iatoms_pairs.size();
        all_cells.reserve(group_size);

        for (int i = 0; i < group_size; i++){
            std::string group = gname_iatoms_pairs[i].first;
            std::vector<int> iatoms = gname_iatoms_pairs[i].second;
            int num_iatoms = iatoms.size();

            // set root cell
            All_cells all_cell;
            all_cell.group = group;
            all_cell.cells = std::vector<Cell>();
            // all_cell.cells.reserve(num_iatoms);

            std::vector<Cell>& cells = all_cell.cells;
            cells.push_back(Cell(n_crit));
            calculate_rc(iatoms);
            cells[0].rc = result_root.rc;
            cells[0].r = result_root.r;

            for (int j = 0; j < num_iatoms; j++){

                int current_cell = 0;
                int index_atom = iatoms[j] - 1;

                while (cells[current_cell].nleaf > n_crit) {
                    // cells[current_cell].leaf.push_back(iatoms[j]);
                    cells[current_cell].nleaf += 1;
                    int octant = (t_crd(index_atom, 0) > cells[current_cell].rc(0)) + \
                                 ((t_crd(index_atom, 1) > cells[current_cell].rc(1)) << 1) + \
                                 ((t_crd(index_atom, 2) > cells[current_cell].rc(2)) << 2);

                                 
                    if (!(cells[current_cell].nchild & (1 << octant))){
                        add_child(octant, current_cell, cells);
                    }

                    current_cell = cells[current_cell].child[octant];
                    // std::cerr << "current cell r: " << cells[current_cell].r << std::endl;
                }

                cells[current_cell].leaf.push_back(iatoms[j]);
                cells[current_cell].nleaf += 1;

                if (cells[current_cell].nleaf > n_crit){
                    split_cell(current_cell, cells);
                }
            }
            all_cells.push_back(all_cell);
            // std::cerr << "root cell r: " << cells[0].r << std::endl;
            // std::cerr << "group: " << group << std::endl;
            // std::cerr << "" << std::endl;
        }
        // std::cerr << "number of all cells: " << all_cells[0].cells.size() << std::endl;
        // for (int i = 0; i < all_cells[0].cells.size(); i++){
        //     std::cerr << "number of all cells: " << all_cells[0].cells[i].leaf << std::endl;
        // }
        // print_all_cells();
        int t1 = time_now();
        std::cerr << "setup_all_cells time: " << t1 - t0 << " seconds" << std::endl;
        return all_cells;
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void get_all_cells(const std::vector<All_cells>& all_cells_input){

        all_cells = all_cells_input;  
    };    

    // calculate multipole
    void cal_multipole(VectorXd& multipole, MatrixXd& multipole_j, Vector3d& rc, std::vector<int>& atoms){
        
        int size_atoms = atoms.size();
        for (int i = 0; i < size_atoms; i++){

            int index_atom = atoms[i] - 1;

            Vector3d crd_atom = t_crd.row(index_atom);
            double dx = rc(0) - crd_atom(0);
            double dy = rc(1) - crd_atom(1);
            double dz = rc(2) - crd_atom(2);
            double qj = charges[index_atom];

            double qjdx = qj * dx;
            double qjdy = qj * dy;
            double qjdz = qj * dz;
            
            VectorXd multipole_atom = VectorXd::Zero(10);

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

            for (int j = 0; j < 3; j++){
                multipole_j.row(j) += multipole_atom.transpose() * crd_atom(j);
            }
            // std::cerr << "atom: " << atoms[i] << std::endl;
            // std::cerr << "atom crd: " << crd_atom.transpose() << std::endl;
            // std::cerr << "atom rc: " << rc.transpose() << std::endl;
            // std::cerr << "atom charge: " << qj << std::endl;
            // std::cerr << "dx: " << dx << ", dy: " << dy << ", dz: " << dz << std::endl;
            // std::cerr << "qjdx: " << qjdx << ", qjdy: " << qjdy << ", qjdz: " << qjdz << std::endl;
            // std::cerr << "multipole: " << multipole.transpose() << std::endl;
            // std::cerr << "   " << std::endl;
        }
        // std::cerr << atoms.transpose() <<  multipole.transpose() << std::endl;
    };

    void cal_p(){
        int size_all_cells = all_cells.size();
        for (int i = 0; i < size_all_cells; i++){
            std::vector<Cell>& cells = all_cells[i].cells;
            int size_cells = cells.size();

            for (int j = 0; j < size_cells; j++){
               
                Cell& cell = cells[j];
                if (cell.nchild != 0){
                    continue;
                }
                else {
                    // std::cerr << "calculating multipole for cell: " << j << std::endl;
                    // std::cerr << "" << std::endl;
                    cal_multipole(cell.multipole, cell.multipole_j, cell.rc, cell.leaf);
                }
            }
            // std::cerr << "multipole(cal_p): " << cells[size_cells-1].multipole.transpose() << std::endl;
        }
        // print_multipoles();
    };

    void cal_M2M(){
        int size_all_cells = all_cells.size();
        for (int i = 0; i < size_all_cells; i++){
            std::vector<Cell>& cells = all_cells[i].cells;
            
            int cells_size = cells.size();
            for (int i = 0; i < cells_size - 1; i++){
                int i_inv = cells_size - 1 - i;

                int c = i_inv;
                int p = cells[i_inv].parent;
            
                VectorXd& p_potential = cells[p].multipole;
                VectorXd& c_potential = cells[c].multipole;
                Vector3d& c_rc = cells[c].rc;
                Vector3d& p_rc = cells[p].rc;
                
                double dx = p_rc(0) - c_rc(0);
                double dy = p_rc(1) - c_rc(1);
                double dz = p_rc(2) - c_rc(2);
                double Mx = c_potential(0) * dx;
                double My = c_potential(0) * dy;
                double Mz = c_potential(0) * dz;

                p_potential(0) += c_potential(0);
                p_potential(1) += c_potential(1) + Mx;
                p_potential(2) += c_potential(2) + My;
                p_potential(3) += c_potential(3) + Mz;
                p_potential(4) += c_potential(4) + dx * c_potential(1) + 0.5 * Mx * dx;
                p_potential(5) += c_potential(5) + dy * c_potential(2) + 0.5 * My * dy;
                p_potential(6) += c_potential(6) + dz * c_potential(3) + 0.5 * Mz * dz;
                p_potential(7) += c_potential(7) + dy * c_potential(1) + dx * c_potential(2) + Mx * dy;
                p_potential(8) += c_potential(8) + dz * c_potential(2) + dy * c_potential(3) + My * dz;
                p_potential(9) += c_potential(9) + dx * c_potential(3) + dz * c_potential(1) + Mz * dx;

                MatrixXd& p_potential_j = cells[p].multipole_j;
                MatrixXd& c_potential_j = cells[c].multipole_j;

                for (int j = 0; j < 3; j++){
                    double Mjx = c_potential_j(j, 0) * dx;
                    double Mjy = c_potential_j(j, 0) * dy;
                    double Mjz = c_potential_j(j, 0) * dz;

                    p_potential_j(j, 0) += c_potential_j(j, 0);
                    p_potential_j(j, 1) += c_potential_j(j, 1) + Mjx;
                    p_potential_j(j, 2) += c_potential_j(j, 2) + Mjy;
                    p_potential_j(j, 3) += c_potential_j(j, 3) + Mjz;
                    p_potential_j(j, 4) += c_potential_j(j, 4) + dx * c_potential_j(j, 1) + 0.5 * Mjx * dx;
                    p_potential_j(j, 5) += c_potential_j(j, 5) + dy * c_potential_j(j, 2) + 0.5 * Mjy * dy;
                    p_potential_j(j, 6) += c_potential_j(j, 6) + dz * c_potential_j(j, 3) + 0.5 * Mjz * dz;
                    p_potential_j(j, 7) += c_potential_j(j, 7) + dy * c_potential_j(j, 1) + dx * c_potential_j(j, 2) + Mjx * dy;
                    p_potential_j(j, 8) += c_potential_j(j, 8) + dz * c_potential_j(j, 2) + dy * c_potential_j(j, 3) + Mjy * dz;
                    p_potential_j(j, 9) += c_potential_j(j, 9) + dx * c_potential_j(j, 3) + dz * c_potential_j(j, 1) + Mjz * dx;
                }

                // std::cerr << "calculatiing child " << i_inv << " to parent " << p << std::endl;
                // std::cerr << "child multipole: " << c_potential.transpose() << std::endl;
                // std::cerr << "p_rc: " << p_rc.transpose() << std::endl;
                // std::cerr << "c_rc: " << c_rc.transpose() << std::endl;
                // std::cerr << "dx: " << dx << ", dy: " << dy << ", dz: " << dz << std::endl;
                // std::cerr << "c_potential 0: " << c_potential(0) << std::endl;
                // std::cerr << "Mx: " << Mx << ", My: " << My << ", Mz: " << Mz << std::endl;
                // std::cerr << "cell " << p << " p multipole: " << p_potential.transpose() << std::endl;
                // std::cerr << "   " << std::endl;
            }
        }
        // std::cerr << "multipole: " << all_cells[0].cells[0].multipole.transpose() << std::endl;
    };

    void cal_fiJ(int num_source, int p, std::vector<Cell>& cells){

        int t0 = time_now();
        int idx_source = num_source - 1;
        Vector3d vi = t_vel.row(idx_source);

        std::queue<int> cell_queue;
        cell_queue.push(p);

        while (!cell_queue.empty()){
            int p = cell_queue.front();
            cell_queue.pop();
        
            if (cells[p].nleaf > n_crit){

                for (int octant = 0; octant < 8; octant++){
                    
                    if (cells[p].nchild & (1 << octant)) {
                        
                        int t1 = time_now();
                        int c = cells[p].child[octant];
                        Vector3d& rc = cells[c].rc;
                        Vector3d crd_source = t_crd.row(idx_source);
                        
                        double dx = crd_source(0) - rc(0);
                        double dy = crd_source(1) - rc(1);
                        double dz = crd_source(2) - rc(2);
                        double r = sqrt(dx * dx + dy * dy + dz * dz);
                        int t2 = time_now();
                        // std::cerr << "num_source: "<< num_source << " cell: " << p << " child: " << c <<  std::endl;
                        // std::cerr << "crd_source - rc cal time: " << t2-t1 << std::endl;

                        if (cells[c].r > theta * r){
                            // std::cerr << "kaiki to " << c << std::endl;
                            cell_queue.push(c);
                        }
                        else{
                            int t3 = time_now();
                            // for (int l = 0; l < cells[c].nleaf; l++){
                            //     int num_target = cells[c].leaf[l];
                                
                            //     if (is_nonbonded_pair(num_source, num_target) == false){
                            //         std::cerr << "ValueError: Bonded pair in fmm" << std::endl;
                            //         std::cerr << "num_source: " << num_source << "num_target: " << num_target << std::endl;
                            //         std::cerr << "Please change the parameter of fmm." << std::endl;
                            //         std::exit(1);
                            //     }
                            // }
                            VectorXd bJx = VectorXd::Zero(10);
                            VectorXd bJy = VectorXd::Zero(10);
                            VectorXd bJz = VectorXd::Zero(10);

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

                            int t4 = time_now();

                            // calculate potential
                            VectorXd& potential = cells[c].multipole;
                            double charge = charges[idx_source];
                            double coeff = 332.05221729;
                            VectorXd pot = potential * charge * coeff;

                            double fx = pot.dot(bJx);
                            double fy = pot.dot(bJy);
                            double fz = pot.dot(bJz);
                            Vector3d f = Vector3d(-fx, -fy, -fz);
                            Vector3d riJ = Vector3d(dx, dy, dz);

                            MatrixXd& potential_j = cells[c].multipole_j;
                            MatrixXd bJ = MatrixXd::Zero(10, 3);
                            bJ.col(0) = bJx;
                            bJ.col(1) = bJy;
                            bJ.col(2) = bJz;
                            Matrix3d pot_j = -potential_j * bJ * coeff * charge;

                            int igrp = iatom_to_igroup(idx_source);
                            int Jgrp = iatom_to_igroup(cells[c].leaf[0] - 1);
                            // VectorXd vJ = VectorXd::Zero(3);
                            // for (int l = 0; l < cells[c].nleaf; l++){
                            //     int num_leaf = cells[c].leaf[l];
                            //     if (num_leaf == num_source){
                            //         continue;
                            //     }
                            //     vJ += t_vel.row(num_leaf - 1);

                            // }
                            // vi = vi + vJ;

                            // std::cerr << "cellwise calculation " << std::endl;
                            // std::cerr << "num_source: " << num_source << std::endl;
                            // std::cerr << "num_target: [";
                            // for (size_t l = 0; l < cells[c].leaf.size(); ++l) {
                            //     std::cerr << cells[c].leaf[l];
                            //     if (l < cells[c].leaf.size() - 1) {
                            //         std::cerr << ", ";
                            //     }
                            // };

                            // std::cerr << "]" << std::endl;
                            // std::cerr << "crd_source: " << crd_source.transpose() << std::endl;
                            // std::cerr << "crd_target: " << rc.transpose() << std::endl;
                            // std::cerr << "riJ: " << riJ.transpose() << std::endl;
                            // std::cerr << "r: " << r << std::endl;
                            // std::cerr << "bJx: " << bJx.transpose() << std::endl;
                            // std::cerr << "bJy: " << bJy.transpose() << std::endl;
                            // std::cerr << "bJz: " << bJz.transpose() << std::endl;
                            // std::cerr << "potential: " << potential.transpose() << std::endl;
                            // std::cerr << "charge: " << charge << std::endl;
                            // std::cerr << "pot: " << pot.transpose() << std::endl;
                            // std::cerr << "f: " << f.transpose() << std::endl;
                            // std::cerr << "vi: " << vi.transpose() << std::endl;
                            // std::cerr << "riJ: " << riJ.transpose() << std::endl;
                            // std::cerr << "igrp: " << igrp << std::endl;
                            // std::cerr << "Jgrp: " << Jgrp << std::endl;
                            // // std::cerr << "vJ: " << vJ.transpose() << std::endl;
                            // std::cerr << "   " << std::endl;
                            cal_flux_cellwise(f, vi, riJ, igrp, Jgrp);

                            int t5 = time_now();
                            // std::cerr << "cal_fiJ prep time: " << t4 - t3 << std::endl;
                            // std::cerr << "cal_fiJ time: " << t5 - t4 << std::endl;
                        }      
                    }
                }
            }
            else {
                int t6 = time_now();
                for (int l = 0; l < cells[p].nleaf; l++){
                    int num_target = cells[p].leaf[l];
                    int idx_target = num_target - 1;

                    if (num_target == num_source){
                        continue;
                    }
                    if (is_nonbonded_pair(num_source, num_target) == false){
                        continue;
                    }

                    Vector3d crd_source = t_crd.row(idx_source);
                    Vector3d crd_target = t_crd.row(idx_target);
                    Vector3d rij = crd_source - crd_target;
                    
                    double r = sqrt(rij.dot(rij));
                    double inv_r = 1.0 / r;
                    double coeff = 332.05221729;

                    double charge_i = charges[idx_source];
                    double charge_j = charges[idx_target];
                    double qij = coeff * charge_i * charge_j;
                    qij = qij * inv_r * inv_r * inv_r;
                
                    Vector3d fij = rij * qij;

                    int igrp = iatom_to_igroup(idx_source);
                    int jgrp = iatom_to_igroup(idx_target);

                    // std::cerr << "atomwise calculation " << std::endl;
                    // std::cerr << "num_source: " << num_source << std::endl;
                    // std::cerr << "num_target: " << num_target << std::endl;
                    // std::cerr << "crd_source: " << crd_source.transpose() << std::endl;
                    // std::cerr << "crd_target: " << crd_target.transpose() << std::endl;
                    // std::cerr << "rij: " << rij.transpose() << std::endl;
                    // std::cerr << "vi: " << vi.transpose() << std::endl;
                    // std::cerr << "vj: " << vj.transpose() << std::endl;
                    // std::cerr << "vij: " << vij.transpose() << std::endl;
                    // std::cerr << "r: " << r << std::endl;
                    // std::cerr << "charge_i: " << charge_i << std::endl;
                    // std::cerr << "charge_j: " << charge_j << std::endl;
                    // std::cerr << "qij: " << qij << std::endl;
                    // std::cerr << "fij: " << fij.transpose() << std::endl;
                    // std::cerr << "" << std::endl;

                    cal_flux_atomwise(fij, vi, rij, igrp, jgrp);

                }
                int t7 = time_now();
                // std::cerr << "cal_fij time(/pair): " << (t7 - t6)/cells[p].nleaf << std::endl;  
            }
        }
    };

    void cal_hflux_cellwise(Vector3d fiJ, Vector3d vi, Vector3d r, Matrix3d pot_j, int igrp, int Jgrp){
        
        Vector3d h_ij = r * (fiJ.dot(vi)) * 0.5 - pot_j * vi * 0.5;

        if (igrp != Jgrp){
            hflux_ij[igrp-1][Jgrp-1] += h_ij;
            hflux_ij[Jgrp-1][igrp-1] += h_ij;
        }
        else{
            hflux_ij[igrp-1][Jgrp-1] += 0.5 * h_ij;
        }
        count_cell += 1;
    };

    void cal_hflux_atomwise(Vector3d fij, Vector3d vi, Vector3d rij, int igrp, int jgrp){
        
        Vector3d h_ij = rij * (fij.dot(vi)) * 0.5;

        if (igrp != jgrp){
            hflux_ij[igrp-1][jgrp-1] += h_ij;
            hflux_ij[jgrp-1][igrp-1] += h_ij;
        }
        else{
            hflux_ij[igrp-1][jgrp-1] += 0.5 * h_ij;
        }
        count_atom += 1;
    };

    void cal_eflux_cellwise(Vector3d fiJ, Vector3d vi, Vector3d r, Matrix3d pot_j, int igrp, int Jgrp){
        
        double e_ij = fiJ.dot(vi) * 0.5;

        if (igrp != Jgrp){
            eflux_ij(igrp-1, Jgrp-1) += e_ij;
            eflux_ij(Jgrp-1, igrp-1) += -e_ij;
        }
        else{
            eflux_ij(igrp-1, Jgrp-1) += 0.5 * e_ij;         // it is not correct
        }
        count_cell += 1;
    };

    void cal_eflux_atomwise(Vector3d fij, Vector3d vi, Vector3d rij, int igrp, int jgrp){

        double e_ij = fij.dot(vi) * 0.5;

        if (igrp != jgrp){
            eflux_ij(igrp-1, jgrp-1) += e_ij;
            eflux_ij(jgrp-1, igrp-1) += -e_ij;
        }
        else{
            eflux_ij(igrp-1, jgrp-1) += e_ij;
        }
        count_atom += 1;
    };

    // check if the pairs are in the interact table
    bool is_nonbonded_pair(int source, int target){

        bool is_nonbond = false;
        if (source > target){
            std::swap(source, target);
        }
        for (const auto& table : interact_table){
            for (const auto& tuple : table){
                int i = std::get<0>(tuple);
                int j_beg = std::get<1>(tuple);
                int j_end = std::get<2>(tuple);
                if (i == source && (j_beg <= target && target <= j_end)){
                    is_nonbond = true;
                    break;
                }
            }
        }
        return is_nonbond;
    };

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

    void cal_force(){

        int size_table = gpair_table.size();
        for (int i = 0; i < size_table; i++){

            std::string source = gpair_table[i].first;
            std::vector<std::string> targets = gpair_table[i].second;
            int size_targets = targets.size();
        
            std::vector<int> source_atoms = get_atoms(source, gname_iatoms_pairs);
            int source_size = source_atoms.size();

            for (int j = 0; j < source_size; j++){

                int num_source = source_atoms[j];

                for (int k = 0; k < size_targets; k++){      

                    std::string target = targets[k];
                    std::vector<Cell> cells = get_cells(target, all_cells);
                    cal_fiJ(num_source, 0, cells);
                }
            }

            std::vector<Cell> cells_source = get_cells(source, all_cells);
            for (int l = 0; l < size_targets; l++){

                std::string target = targets[l];
                std::vector<int> target_atoms = get_atoms(target, gname_iatoms_pairs);
                int target_size = target_atoms.size();

                for (int m = 0; m < target_size; m++){
                    int num_target = target_atoms[m];
                    cal_fiJ(num_target, 0, cells_source);
                }
            }
        }
    };

    void cal_coulomb_flux_fmm(const std::vector<All_cells>& all_cells) {
        int t0 = time_now();
        get_all_cells(all_cells);
        int t1 = time_now();
        cal_p();
        int t2 = time_now();
        cal_M2M();
        int t3 = time_now();
        cal_force();
        int t4 = time_now();
        std::cerr << "cal_p time: " << t2 - t1 << std::endl;
        std::cerr << "cal_M2M time: " << t3 - t2 << std::endl;
        std::cerr << "cal_force time: " << t4 - t3 << std::endl;
        std::cerr << "total time(cal_p -> cal_force): " << t4 - t1 << std::endl;
        std::cerr << "   " << std::endl;
    };

    int time_now(){
 
        auto start = std::chrono::system_clock::now();
        // Simulate some work
        auto now = std::chrono::duration_cast<std::chrono::microseconds>(start.time_since_epoch()).count();
        // float nnow = now * 0.000001;
        return now;
    };

    // for debug
    void check_cells(){

        std::cerr << "start check cells" << std::endl;

        int size_all_cells = all_cells.size();
        for (int i = 0; i < size_all_cells; i++){
            std::cerr << "group: " << all_cells[i].group << std::endl;
            std::vector<int> atoms = get_atoms(all_cells[i].group, gname_iatoms_pairs);
            int size_atoms = atoms.size();

            std::vector<cal_fmm::Cell> cells = all_cells[i].cells;
            int size_cells = cells.size();

            // 1. check if atoms are in the cell
            for (int j = 0; j < size_cells; j++){
                
                Vector3d rc = cells[j].rc;                      // center of the cell
                double r = cells[j].r;                          // radius of the cell  
                Vector3d ra = Vector3d(r, r, r);
                Vector3d rmin = rc - ra;                        // minimum coordinate of the cell
                Vector3d rmax = rc + ra;                        // maximum coordinate of the cell

                for (int k = 0; k < cells[j].nleaf; k++){

                    int atom = cells[j].leaf[k];
                    if (atom <= 0){
                        std::cerr << "atom=0 " << atom << std::endl;
                        continue;
                    }
                    Vector3d crd_atom = t_crd.row(atom-1);

                    // check if the atom is in the cell
                    if (crd_atom(0) < rmin(0) || crd_atom(0) > rmax(0) || \
                        crd_atom(1) < rmin(1) || crd_atom(1) > rmax(1) || \
                        crd_atom(2) < rmin(2) || crd_atom(2) > rmax(2)){

                        
                        std::cerr << "error: " << atom << " " << crd_atom.transpose() << std::endl;
                        std::cerr << "cell_index: " << j << std::endl;
                        std::cerr << "rmin: " << rmin.transpose() << std::endl;
                        std::cerr << "rmax: " << rmax.transpose() << std::endl;
                        if (crd_atom(0) < rmin(0) || crd_atom(0) > rmax(0)){
                            std::cerr << "atom x out of range" << std::endl;
                        }
                        if (crd_atom(1) < rmin(1) || crd_atom(1) > rmax(1)){
                            std::cerr << "atom y out of range" << std::endl;
                        }
                        if (crd_atom(2) < rmin(2) || crd_atom(2) > rmax(2)){
                            std::cerr << "atom z out of range" << std::endl;
                        }
                        std::cerr << "    " << std::endl;
                    }
                }
            }

            for (int j = 0; j < size_cells; j++){

                std::vector<int> atoms_inc = std::vector<int>();
                Vector3d rc = cells[j].rc;                      // center of the cell
                double r = cells[j].r;                          // radius of the cell  
                Vector3d ra = Vector3d(r, r, r);
                Vector3d rmin = rc - ra;                        // minimum coordinate of the cell
                Vector3d rmax = rc + ra;                        // maximum coordinate of the cell

                // check if the cell is split correctly
                for (int k = 0; k < size_atoms; k++){
                    int atom = atoms[k];
                    if (atom <= 0){
                        std::cerr << "atom=0 " << atom << std::endl;
                        continue;
                    }
                    Vector3d crd_atom = t_crd.row(atom-1);
                    if (crd_atom(0) >= rmin(0) && crd_atom(0) <= rmax(0) && \
                        crd_atom(1) >= rmin(1) && crd_atom(1) <= rmax(1) && \
                        crd_atom(2) >= rmin(2) && crd_atom(2) <= rmax(2)){
                        atoms_inc.push_back(atom);
                    }
                }

                // check if the number of atoms in the cell is correct
                int atoms_inc_size = atoms_inc.size();
                if (cells[j].nleaf != atoms_inc_size){
                    std::cerr << "error: cell " << j << std::endl;
                    std::cerr << "num of leaves: " << cells[j].nleaf << std::endl;
                    std::cerr << "atoms_included: " << atoms_inc.size() << std::endl;
                }

                // check if the atoms in the cell are correct
                if (cells[j].nchild == 0){
                    for (int k = 0; k < cells[j].nleaf; k++){

                        int atom_k = cells[j].leaf[k];
                        if (atom_k != atoms_inc[k]){
                            std::cerr << "error: cell " << j << std::endl;
                            std::cerr << "atom(in all_cells): " << atom_k << std::endl;
                            std::cerr << "atoms_inc: " << atoms_inc[k] << std::endl;
                            std::cerr << "atom_k: [";

                            for (size_t l = 0; l < cells[j].leaf.size(); ++l) {
                                std::cerr << cells[j].leaf[l];
                                if (l < cells[j].leaf.size() - 1) {
                                    std::cerr << ", ";
                                }
                            };
                            std::cerr << "]" << std::endl;

                            std::cerr << "atoms_inc: " ;
                            for (const auto& atom : atoms_inc) {
                                std::cerr << atom << " ";
                            }
                            std::cerr << std::endl;
                            std::cerr << "  " << std::endl;
                        }
                    }
                }
                else {
                    for (int k = 0; k < n_crit; k++){
                        int atom_k = cells[j].leaf[k];
                        if (atom_k != atoms_inc[k]){
                            std::cerr << "error: cell " << j << std::endl;
                            std::cerr << "atom(in all_cells): " << atom_k << std::endl;
                            std::cerr << "atoms_inc: " << atoms_inc[k] << std::endl;
                            std::cerr << "atom_k: [";

                            for (size_t l = 0; l < cells[j].leaf.size(); ++l) {
                                std::cerr << cells[j].leaf[l];
                                if (l < cells[j].leaf.size() - 1) {
                                    std::cerr << ", ";
                                }
                            };
                            std::cerr << "]" << std::endl;

                            std::cerr << "atoms_inc: " ;
                            for (const auto& atom : atoms_inc) {
                                std::cerr << atom << " ";
                            }
                            std::cerr << std::endl;
                            std::cerr << "  " << std::endl;
                        }
                    }
                }

                // check if the parent and child are correct
                if (cells[j].nchild != 0){
                    for (int k = 0; k < 8; k++){
                        if (cells[j].child[k] == 0){
                            continue;
                        }
                        int child = cells[j].child[k];
                        if (cells[child].parent != j){
                            std::cerr << "error: cell " << j << std::endl;
                            std::cerr << "child: " << child << std::endl;
                            std::cerr << "parent: " << cells[child].parent << std::endl;
                        }

                        // check if the radius of the child is correct
                        if (cells[j].r != cells[child].r * 2.0){
                            std::cerr << "error: cell " << j << std::endl;
                            std::cerr << "child: " << child << std::endl;
                            std::cerr << "r: " << cells[j].r << std::endl;
                            std::cerr << "r_child: " << cells[child].r * 2.0 << std::endl;
                        }
                        
                    }
                }
            }

            // check if M2M is correct
            for (int j = 0; j < size_cells; j++){

                Cell cell = cells[j];
                VectorXd multi_test = VectorXd::Zero(10);
                Vector3d rcenter = cell.rc;
                for (int k = 0; k < cell.nleaf; k++){

                    int index_leaf = cell.leaf[k] - 1;
                    Vector3d crd_leaf = t_crd.row(index_leaf);
                    Vector3d r = rcenter - crd_leaf;
                    double charge = charges[index_leaf];

                    multi_test(0) += charge * 1.0;
                    multi_test(1) += charge * r(0);
                    multi_test(2) += charge * r(1);
                    multi_test(3) += charge * r(2);
                    multi_test(4) += charge * r(0) * r(0) * 0.5;
                    multi_test(5) += charge * r(1) * r(1) * 0.5;
                    multi_test(6) += charge * r(2) * r(2) * 0.5;
                    multi_test(7) += charge * r(0) * r(1);
                    multi_test(8) += charge * r(1) * r(2);
                    multi_test(9) += charge * r(2) * r(0);
                } 
                VectorXd& multipole = cell.multipole;
                VectorXd diff = multi_test - multipole;
                int size_diff = diff.size();
                VectorXd norm_diff = VectorXd::Zero(10);
                for (int k = 0; k < diff.size(); k++){
                    norm_diff(k) = diff(k) / multipole(k);
                }
                if (multi_test != multipole){
                    std::cerr << "M2M is not equal in cell " << j << std::endl;
                    std::cerr << "multipole: " << multipole.transpose() << std::endl;
                    std::cerr << "multi_test: " << multi_test.transpose() << std::endl;
                    std::cerr << "diff: " << diff.transpose() << std::endl;
                    std::cerr << "normalized diff: " << norm_diff.transpose() << std::endl;
                }
            }
        }
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
        .def_readwrite("interact_table", &cal_fmm::interact_table)
        .def_readwrite("gpair_table", &cal_fmm::gpair_table)
        .def_readwrite("gname_iatoms_pairs", &cal_fmm::gname_iatoms_pairs)
        .def_readwrite("iatom_to_igroup", &cal_fmm::iatom_to_igroup)
        .def_readwrite("all_cells", &cal_fmm::all_cells)
        .def_readwrite("hflux_ij", &cal_fmm::hflux_ij)
        .def_readwrite("eflux_ij", &cal_fmm::eflux_ij)
        ;
    
    
    py::class_<cal_fmm::Cell>(m, "Cell")
        .def(py::init<int>())
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
}