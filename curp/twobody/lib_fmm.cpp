#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
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
    int n_crit;
    float theta;
    std::vector<double> charges;
    MatrixXd t_crd;

    std::vector<std::pair<std::string, std::vector<std::string>>> gpair_table_fmm;
    std::vector<std::pair<std::string, std::vector<int>>> gname_iatoms_pairs;

    std::vector<int> atomwise_i;
    std::vector<int> atomwise_j;
    std::vector<Vector3d> atomwise_f;
    std::vector<Vector3d> atomwise_r;

    std::vector<int> cellwise_i;
    std::vector<VectorXi> cellwise_J;
    std::vector<Vector3d> cellwise_f;
    std::vector<Vector3d> cellwise_r;


    // constructor
    cal_fmm():
        natom(0),
        n_crit(0),
        theta(0.0),
        charges(std::vector<double>()),
        t_crd(MatrixXd::Zero(natom, 3)),
        gpair_table_fmm(std::vector<std::pair<std::string, std::vector<std::string>>>()),
        gname_iatoms_pairs(std::vector<std::pair<std::string, std::vector<int>>>())
    {};


    void setup(const int input_natom, const int& input_n_crit, const double& input_theta, const std::vector<double>& input_charges, \
        const std::vector<std::pair<std::string, std::vector<int>>>& input_gname_iatoms_pairs, \
        const std::vector<std::pair<std::string, std::vector<std::string>>>& input_gpair_table_fmm){
        natom  = input_natom;
        n_crit = input_n_crit;
        theta  = input_theta;
        charges = input_charges;
        gname_iatoms_pairs = input_gname_iatoms_pairs;
        gpair_table_fmm = input_gpair_table_fmm;
    };

    // read trajectory
    void initialize(const MatrixXd& crd){
        t_crd = crd;
        atomwise_i = std::vector<int>();
        atomwise_j = std::vector<int>();
        atomwise_f = std::vector<Vector3d>();
        atomwise_r = std::vector<Vector3d>();
        cellwise_i = std::vector<int>();
        cellwise_J = std::vector<VectorXi>();
        cellwise_f = std::vector<Vector3d>();
        cellwise_r = std::vector<Vector3d>();

    };

////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct Cell {
        int       nleaf;
        VectorXi  leaf;
        int       nchild;
        VectorXi  child;
        int       parent;
        Vector3d  rc;
        double    r;
        VectorXd  multipole;

        Cell(int n_crit):
            nleaf(0),                               // number of atoms(leaf) in the cell
            leaf(VectorXi::Zero(n_crit)),           // index of atoms in the cell
            nchild(0),                              // number of child cells
            child(VectorXi::Zero(8)),               // index of 8 child cells
            parent(0),                              // index of parent cell
            rc(Vector3d::Zero()),                   // center of the cell
            r(0.0),                                 // radius of the cell
            multipole(VectorXd::Zero(10))           // 10 multipoles
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
        int size_atoms = distance(atoms.begin(), atoms.end());
        
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

        cells[new_child].r = cells[current_cell].r * 0.5;
        cells[new_child].rc(0) = cells[current_cell].rc(0) + cells[new_child].r * ((octant & 1) * 2 - 1);
        cells[new_child].rc(1) = cells[current_cell].rc(1) + cells[new_child].r * ((octant & 2) - 1);
        cells[new_child].rc(2) = cells[current_cell].rc(2) + cells[new_child].r * ((octant & 4) / 2 - 1);

        cells[new_child].parent = current_cell;
        cells[current_cell].child(octant) = new_child;
        cells[current_cell].nchild = (cells[current_cell].nchild | (1 << octant));
    };

    void split_cell(int current_cell, std::vector<Cell>& cells){

        for (int i = 0; i < cells[current_cell].nleaf; i++){

            int atom_current = cells[current_cell].leaf(i);
            int index_atom_current = atom_current - 1;
            int octant = (t_crd(index_atom_current, 0) > cells[current_cell].rc(0)) + \
                            ((t_crd(index_atom_current, 1) > cells[current_cell].rc(1)) << 1) + \
                            ((t_crd(index_atom_current, 2) > cells[current_cell].rc(2)) << 2);

            if  (!(cells[current_cell].nchild & (1 << octant))){
                add_child(octant, current_cell, cells);

            }
            int child_cell = cells[current_cell].child(octant);
            cells[child_cell].leaf(cells[child_cell].nleaf) = atom_current;
            cells[child_cell].nleaf += 1;

            if (cells[child_cell].nleaf >= n_crit){
                split_cell(child_cell, cells);
            }
        }
    };
        
    std::vector<All_cells> setup_all_cells(){

        all_cells = std::vector<All_cells>();
        int group_size = gname_iatoms_pairs.size();
        for (int i = 0; i < group_size; i++){
            std::string group = gname_iatoms_pairs[i].first;
            std::vector<int> iatoms = gname_iatoms_pairs[i].second;

            // set root cell
            All_cells all_cell;
            all_cell.group = group;
            all_cell.cells = std::vector<Cell>();
            all_cell.cells.push_back(Cell(n_crit));
            calculate_rc(iatoms);
            all_cell.cells[0].rc = result_root.rc;
            all_cell.cells[0].r = result_root.r;

            int num_iatoms = iatoms.size();

            for (int j = 0; j < num_iatoms; j++){

                int current_cell = 0;
                int index_atom = iatoms[j] - 1;

                while (all_cell.cells[current_cell].nleaf >= n_crit) {
                    all_cell.cells[current_cell].nleaf += 1;
                    int octant = (t_crd(index_atom, 0) > all_cell.cells[current_cell].rc(0)) + \
                                 ((t_crd(index_atom, 1) > all_cell.cells[current_cell].rc(1)) << 1) + \
                                 ((t_crd(index_atom, 2) > all_cell.cells[current_cell].rc(2)) << 2);

                                 
                    if (!(all_cell.cells[current_cell].nchild & (1 << octant))){
                        add_child(octant, current_cell, all_cell.cells);
                    }

                    current_cell = all_cell.cells[current_cell].child(octant);

                }
            
                all_cell.cells[current_cell].leaf(all_cell.cells[current_cell].nleaf) = iatoms[j];
                all_cell.cells[current_cell].nleaf += 1;

                if (all_cell.cells[current_cell].nleaf >= n_crit){
                    split_cell(current_cell, all_cell.cells);
                }
            }
            all_cells.push_back(all_cell);
        }
        // std::cerr << "number of all cells: " << all_cells[0].cells.size() << std::endl;
        // for (int i = 0; i < all_cells[0].cells.size(); i++){
        //     std::cerr << "number of all cells: " << all_cells[0].cells[i].leaf.transpose() << std::endl;
        // }
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
                std::cerr << "cell: " << j << " nleaf: " << cells[j].nleaf << " leaves: " << cells[j].leaf.transpose() << std::endl;
                std::cerr << "nchild: " << cells[j].nchild << " child: " << cells[j].child.transpose() << std::endl;
                std::cerr << "parent: " << cells[j].parent << std::endl;
                std::cerr << "rc: " << cells[j].rc.transpose() << std::endl;
                std::cerr << "r: " << cells[j].r << std::endl;
                if (cells[j].multipole(0) != 0){
                    std::cerr << "multipole: " << cells[j].multipole.transpose() << std::endl;
                }
                if (cells[j].nchild != 0 && cells[j].multipole[0] != 0){
                    std::cerr << "error: multipole exists: " << cells[j].child.transpose() << std::endl;
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
    VectorXd cal_multipole(VectorXd& multipole, Vector3d& rc, VectorXi& atoms){
        
        int size_atoms = atoms.size();
        for (int i = 0; i < size_atoms; i++){

            int index_atom = atoms(i) - 1;

            Vector3d crd_atom = t_crd.row(index_atom);
            double dx = rc(0) - crd_atom(0);
            double dy = rc(1) - crd_atom(1);
            double dz = rc(2) - crd_atom(2);
            double qj = charges[index_atom];

            double qjdx = qj * dx;
            double qjdy = qj * dy;
            double qjdz = qj * dz;
            
            multipole(0) = multipole(0) + qj * 1.0;
            multipole(1) = multipole(1) + qjdx;
            multipole(2) = multipole(2) + qjdy;
            multipole(3) = multipole(3) + qjdz;
            multipole(4) = multipole(4) + qjdx * dx;
            multipole(5) = multipole(5) + qjdy * dy;
            multipole(6) = multipole(6) + qjdz * dz;
            multipole(7) = multipole(7) + qjdx * 2 * dy;
            multipole(8) = multipole(8) + qjdy * 2 * dz;
            multipole(9) = multipole(9) + qjdz * 2 * dx;
            
        }
        return multipole;
    };

    void cal_p(){
        int size_all_cells = all_cells.size();
        for (int i = 0; i < size_all_cells; i++){
            std::vector<Cell>& cells = all_cells[i].cells;
            int size_cells = cells.size();

            for (int j = 0; j < size_cells; j++){

                if (cells[j].nchild != 0){
                    continue;
                }
                else {
                    cells[j].multipole += cal_multipole(cells[j].multipole, cells[j].rc, cells[j].leaf);
                
                }
            }
            // std::cerr << "multipole: " << cells[size_cells-1].multipole.transpose() << std::endl;
        }
    };

    void cal_M2M(){
        int size_all_cells = all_cells.size();
        for (int i = 0; i < size_all_cells; i++){
            std::vector<Cell>& cells = all_cells[i].cells;
            
            int cells_size = cells.size();
            for (int i = 0; i < cells_size; i++){
                int i_inv = cells_size - 1 - i;

                int c = i_inv;
                int p = cells[i_inv].parent;
            
                VectorXd p_potential = cells[p].multipole;
                VectorXd c_potential = cells[c].multipole;
                VectorXd c_rc = cells[c].rc;
                VectorXd p_rc = cells[p].rc;
                
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

                // std::cerr << "p multipole: " << p_potential.transpose() << std::endl;
                // std::cerr << "c multipole: " << c_potential.transpose() << std::endl;
            }
        }
        // std::cerr << "multipole: " << all_cells[0].cells[0].multipole.transpose() << std::endl;
    };


    void cal_fiJ(int num_source, int p, std::vector<Cell>& cells){

        int idx_source = num_source - 1;

        if (cells[p].nleaf > n_crit){

            for (int octant = 0; octant < 8; octant++){
                
                if (cells[p].nchild & (1 << octant)) {
                    
                    int c = cells[p].child(octant);
                    Vector3d rc = cells[c].rc;
                    Vector3d crd_source = t_crd.row(idx_source);
                    
                    double dx = crd_source(0) - rc(0);
                    double dy = crd_source(1) - rc(1);
                    double dz = crd_source(2) - rc(2);
                    double r = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));

                    if (cells[c].r> theta * r){
                        cal_fiJ(num_source, c, cells);
                    }
                    else{

                        VectorXd bJx = VectorXd::Zero(10);
                        VectorXd bJy = VectorXd::Zero(10);
                        VectorXd bJz = VectorXd::Zero(10);

                        double inv_r = 1.0 / r;
                        double r2 = inv_r * inv_r;
                        double r3 = r2 * inv_r;
                        double r5 = r3 * r2;
                        double r7 = r5 * r2;

                        double dx2 = dx * dx;
                        double dy2 = dy * dy;
                        double dz2 = dz * dz;

                        double dxdy = dx * dy;
                        double dydz = dy * dz;
                        double dzdx = dz * dx;

                        double dxr5 = 3 * dx * r5;
                        double dyr5 = 3 * dy * r5;
                        double dzr5 = 3 * dz * r5;

                        double dxdydz = 15 * dxdy * dz * r7;
                        double dx2dy  = 15 * dx2 * dy * r7;
                        double dy2dz  = 15 * dy2 * dz * r7;
                        double dz2dx  = 15 * dz2 * dx * r7;
                        double dy2dx  = 15 * dy2 * dx * r7;
                        double dz2dy  = 15 * dz2 * dy * r7;
                        double dx2dz  = 15 * dx2 * dz * r7;


                        // calculate bJx
                        bJx(0) = -dx * r3;                          // -dx/r^3
                        bJx(1) = -r3 + 3 * dx2 * r5;                // -1/r^3  + 3dx^2/r^5
                        bJx(2) = 3 * dxdy * r5;                     // 0       + 3dydx/r^5
                        bJx(3) = 3 * dzdx * r5;                     // 0       + 3dzdx/r^5
                        bJx(4) = 3 * dxr5 - 15 * dx2 * dx * r7;     // 9dx/r^5 - 15dx^3/r^7
                        bJx(5) =     dxr5 - dy2dx;                  // 3dx/r^5 - 15dy^2dx/r^7
                        bJx(6) =     dxr5 - dz2dx;                  // 3dx/r^5 - 15dz^2dx/r^7
                        bJx(7) =     dyr5 - dx2dy;                  // 3dy/r^5 - 15dxdy/r^7
                        bJx(8) =          - dxdydz;                 // 0       - 15dydzdx/r^7
                        bJx(9) =     dzr5 - dx2dz;                  // 3dz/r^5 - 15dx^2dz/r^7

                        // calculate bJy
                        bJy(0) = -dy * r3;                          // -dy/r^3
                        bJy(1) = bJx(2);                            // 0       + 3dxdy/r^5
                        bJy(2) = -r3 + 3 * dy2 * r5;                // -1/r^3  + 3dy^2/r^5
                        bJy(3) =     dydz * r5;                     // 0       + 3dzdy/r^5
                        bJy(4) =     dyr5 - dx2dy;                  // 3dy/r^5 - 15dx^2dy/r^7
                        bJy(5) = 3 * dyr5 - 15 * dy2 * dy * r7;     // 9dy/r^5 - 15dy^3/r^7
                        bJy(6) =     dyr5 - dz2dy;                  // 3dy/r^5 - 15dz^2dy/r^7
                        bJy(7) =     dxr5 - dy2dx;                  // 3dx/r^5 - 15dxdy^2/r^7
                        bJy(8) =     dzr5 - dy2dz;                  // 3dz/r^5 - 15dydz^2/r^7
                        bJy(9) =         - dxdydz;                  // 0       - 15dzdxdy/r^7

                        // calculate bJz
                        bJz(0) = -dz * r3;                          // -dz/r^3
                        bJz(1) = bJx(3);                            // 0       + 3dxdz/r^5
                        bJz(2) = bJy(3);                            // 0       + 3dydz/r^5
                        bJz(3) = -r3 + 3 * dz2 * r5;                // -1/r^3  + 3dz^2/r^5
                        bJz(4) =     dzr5 - dx2dz;                  // 3dz/r^5 - 15dx^2dz/r^7
                        bJz(5) =     dzr5 - dy2dz;                  // 3dz/r^5 - 15dy^2dz/r^7
                        bJz(6) = 3 * dzr5 - 15 * dz2 * dz * r7;     // 9dz/r^5 - 15dz^3/r^7
                        bJz(7) =          - dxdydz;                 // 0       - 15dxdydz/r^7
                        bJz(8) =     dyr5 - dz2dy;                  // 3dy/r^5 - 15dydz^2/r^7
                        bJz(9) =     dxr5 - dz2dx;                  // 3dx/r^5 - 15dzdxdz/r^7

                        // calculate potential
                        VectorXd potential = cells[c].multipole;
                        double charge = charges[idx_source];
                        potential = potential * charge;

                        double fx = potential.dot(bJx);
                        double fy = potential.dot(bJy);
                        double fz = potential.dot(bJz);
                        Vector3d f = Vector3d(fx, fy, fz);
                        Vector3d rij = Vector3d(dx, dy, dz);

                        cellwise_i.push_back(num_source);
                        cellwise_J.push_back(cells[c].leaf);
                        cellwise_f.push_back(f);
                        cellwise_r.push_back(rij);
                        // std::cerr << "atomwise_i: " << num_source << " atomwise_j: " << cellwise_J[0] << std::endl;
                        // std::cerr << "atomwise_f: " << f.transpose() << std::endl;
        
                    }      
                }
            }
        }
        else {
            for (int l; l < cells[p].nleaf; l++){
                int num_target = cells[p].leaf(l);
                int idx_target = num_target - 1;

                if (num_target == num_source){
                    continue;
                }

                Vector3d crd_source = t_crd.row(idx_source);
                Vector3d crd_target = t_crd.row(idx_target);

                Vector3d rij = crd_source - crd_target;
                double r = sqrt(pow(rij(0), 2) + pow(rij(1), 2) + pow(rij(2), 2));
                double inv_r = 1.0 / r;
                double coeff = 332.05221729;

                double charge_i = charges[idx_source];
                double charge_j = charges[idx_target];
                double qij = coeff * charge_i * charge_j;
                qij = qij * inv_r * inv_r * inv_r;
            
                Vector3d Fij = rij * qij;

                atomwise_i.push_back(num_source);
                atomwise_j.push_back(num_target);
                atomwise_f.push_back(Fij);
                atomwise_r.push_back(rij);
            }
        }
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

        int size_table = gpair_table_fmm.size();
        for (int i = 0; i < size_table; i++){

            std::string source = gpair_table_fmm[i].first;
            std::vector<std::string> targets = gpair_table_fmm[i].second;
        
            std::vector<int> source_atoms = get_atoms(source, gname_iatoms_pairs);
            int source_size = source_atoms.size();

            for (int j = 0; j < source_size; j++){

                int num_source = source_atoms[j];
                int size_targets = targets.size();

                for (int k = 0; k < size_targets; k++){      

                    std::string target = targets[k];
                    std::vector<Cell> cells = get_cells(target, all_cells);
                    cal_fiJ(num_source, 0, cells);
                }
            }
        }
    };

    void cal_force_fmm(){
        cal_p();
        cal_M2M();
        cal_force();
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

                if (cells[j].nchild != 0){
                    continue;
                }
                else {
                    Vector3d rc = cells[j].rc;                      // center of the cell
                    double r = cells[j].r;                          // radius of the cell  
                    Vector3d ra = Vector3d(r, r, r);
                    Vector3d rmin = rc - ra;                        // minimum coordinate of the cell
                    Vector3d rmax = rc + ra;                        // maximum coordinate of the cell

                    for (int k = 0; k < cells[j].nleaf; k++){

                        int atom = cells[j].leaf(k);
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
                if (cells[j].nleaf != atoms_inc.size()){
                    std::cerr << "error: cell " << j << std::endl;
                    std::cerr << "num of leaves: " << cells[j].nleaf << std::endl;
                    std::cerr << "atoms_included: " << atoms_inc.size() << std::endl;
                }

                // check if the atoms in the cell are correct
                if (cells[j].nchild == 0){
                    for (int k = 0; k < cells[j].nleaf; k++){

                        int atom_k = cells[j].leaf(k);
                        if (atom_k != atoms_inc[k]){
                            std::cerr << "error: cell " << j << std::endl;
                            std::cerr << "atom(in all_cells): " << atom_k << std::endl;
                            std::cerr << "atoms_inc: " << atoms_inc[k] << std::endl;
                            std::cerr << "atom_k: " << cells[j].leaf.transpose() << std::endl;
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
                        int atom_k = cells[j].leaf(k);
                        if (atom_k != atoms_inc[k]){
                            std::cerr << "error: cell " << j << std::endl;
                            std::cerr << "atom(in all_cells): " << atom_k << std::endl;
                            std::cerr << "atoms_inc: " << atoms_inc[k] << std::endl;
                            std::cerr << "atom_k: " << cells[j].leaf.transpose() << std::endl;
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
                        if (cells[j].child(k) == 0){
                            continue;
                        }
                        int child = cells[j].child(k);
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
        }
    };
};

PYBIND11_MODULE(lib_fmm, m){
    
    py::class_<cal_fmm>(m, "cal_fmm")
        .def(py::init<>())
        .def("setup", &cal_fmm::setup)
        .def("initialize", &cal_fmm::initialize)
        .def("calculate_rc", &cal_fmm::calculate_rc)
        .def("setup_all_cells", &cal_fmm::setup_all_cells)
        .def("get_all_cells", &cal_fmm::get_all_cells)
        .def("cal_force_fmm", &cal_fmm::cal_force_fmm)
        .def("cal_p", &cal_fmm::cal_p)
        .def("cal_M2M", &cal_fmm::cal_M2M)
        .def("cal_force", &cal_fmm::cal_force)
        .def("check_cells", &cal_fmm::check_cells)
        .def_readwrite("natom", &cal_fmm::natom)
        .def_readwrite("n_crit", &cal_fmm::n_crit)
        .def_readwrite("theta", &cal_fmm::theta)
        .def_readwrite("charges", &cal_fmm::charges)
        .def_readwrite("t_crd", &cal_fmm::t_crd)
        .def_readwrite("gpair_table_fmm", &cal_fmm::gpair_table_fmm)
        .def_readwrite("gname_iatoms_pairs", &cal_fmm::gname_iatoms_pairs)
        .def_readwrite("atomwise_i", &cal_fmm::atomwise_i)
        .def_readwrite("atomwise_j", &cal_fmm::atomwise_j)
        .def_readwrite("atomwise_f", &cal_fmm::atomwise_f)
        .def_readwrite("atomwise_r", &cal_fmm::atomwise_r)
        .def_readwrite("cellwise_i", &cal_fmm::cellwise_i)
        .def_readwrite("cellwise_J", &cal_fmm::cellwise_J)
        .def_readwrite("cellwise_f", &cal_fmm::cellwise_f)
        .def_readwrite("cellwise_r", &cal_fmm::cellwise_r)
        .def_readwrite("result_root", &cal_fmm::result_root)
        .def_readwrite("all_cells", &cal_fmm::all_cells)
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

    py::class_<cal_fmm::rc_max_r>(m, "rc_max_r")
        .def(py::init<>())
        .def_readwrite("rc", &cal_fmm::rc_max_r::rc)
        .def_readwrite("r", &cal_fmm::rc_max_r::r)
        ;
}