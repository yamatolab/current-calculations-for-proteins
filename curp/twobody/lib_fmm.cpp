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
    // private variables
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


    // create instance
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
            r(0.0),
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
            r(i) = abs(crd_atom.col(i).maxCoeff() - min_r);
            rc(i) = min_r + 0.5 * r(i);
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

    void split_cell(int index_atom, int current_cell, std::vector<Cell>& cells){

        for (int i = 0; i < cells[current_cell].nleaf; i++){

            int index_atom_current = cells[current_cell].leaf(i);
            int octant = (t_crd(index_atom_current, 0) > cells[current_cell].rc(0)) + \
                            ((t_crd(index_atom_current, 1) > cells[current_cell].rc(1)) << 1) + \
                            ((t_crd(index_atom_current, 2) > cells[current_cell].rc(2)) << 2);

            if  (!(cells[current_cell].nchild & (1 << octant))){
                add_child(octant, current_cell, cells);

            }
            int child_cell = cells[current_cell].child(octant);
            cells[child_cell].leaf(cells[child_cell].nleaf) = index_atom_current;
            cells[child_cell].nleaf += 1;
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
                    split_cell(index_atom, current_cell, all_cell.cells);
                }
            }
            all_cells.push_back(all_cell);
        }
        std::cerr << "number of all cells: " << all_cells[0].cells.size() << std::endl;
        for (int i = 0; i < all_cells[0].cells.size(); i++){
            std::cerr << "number of all cells: " << all_cells[0].cells[i].leaf.transpose() << std::endl;
        }
        return all_cells;
    };

////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void get_all_cells(std::vector<All_cells>& all_cells_input){

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
            multipole(9) = multipole(9) + qjdx * 2 * dx;
            
        }
        return multipole;
    };

    void cal_p(int parent_cell){
        int size_all_cells = all_cells.size();
        for (int i = 0; i < size_all_cells; i++){
            std::vector<Cell> cells = all_cells[i].cells;

            if (cells[parent_cell].nleaf >= n_crit){

                for (int j = 0; j < 8; j++){
                    int child_cell = cells[j].nchild;
                    if ((child_cell) & (1 << j)){
                        cal_p(cells[parent_cell].child(child_cell));
                    }
                }
            }
            else {
                cells[parent_cell].multipole += cal_multipole(cells[parent_cell].multipole, cells[parent_cell].rc, cells[parent_cell].leaf);
                std::cerr << "multipole: " << cells[parent_cell].multipole.transpose() << std::endl;
            }
        }
    }

    void cal_M2M(){
        int size_all_cells = all_cells.size();
        for (int i = 0; i < size_all_cells; i++){
            std::vector<Cell> cells = all_cells[i].cells;
            
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

                p_potential(0) = p_potential(0) + c_potential(0);
                p_potential(1) = p_potential(1) + c_potential(0) * dx;
                p_potential(2) = p_potential(2) + c_potential(0) * dy;
                p_potential(3) = p_potential(3) + c_potential(0) * dz;
                p_potential(4) = p_potential(4) + c_potential(1) * dx + 0.5 * c_potential(1) * dx * dx;
                p_potential(5) = p_potential(5) + c_potential(2) * dy + 0.5 * c_potential(2) * dy * dy;
                p_potential(6) = p_potential(6) + c_potential(3) * dz + 0.5 * c_potential(3) * dz * dz;
                p_potential(7) = p_potential(7) + 0.5 * c_potential(2) * dx + 0.5 * c_potential(1) * dx + 0.5 * c_potential(0) * dx * dy;
                p_potential(8) = p_potential(8) + 0.5 * c_potential(3) * dy + 0.5 * c_potential(2) * dy + 0.5 * c_potential(0) * dy * dz;
                p_potential(9) = p_potential(9) + 0.5 * c_potential(1) * dz + 0.5 * c_potential(3) * dz + 0.5 * c_potential(0) * dz * dx;
            }
        }
    };


    void cal_fiJ(int num_source, int p, std::vector<Cell> cells, \
            int idx_cell, int idx_atom){

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
                        cal_fiJ(num_source, c, cells, idx_cell, idx_atom);
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

                        int idx = idx_cell;
                        cellwise[idx].atom_i = num_source;
                        cellwise[idx].atoms_J = cells[c].leaf;
                        cellwise[idx].f = Vector3d(fx, fy, fz);
                        cellwise[idx].r = Vector3d(dx, dy, dz);

                        idx_cell = idx_cell + 1;
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

                int idx = idx_atom;
                atomwise[idx].i = num_source;
                atomwise[idx].j = num_target;
                atomwise[idx].f = Fij;
                atomwise[idx].r = rij;

                idx_atom = idx_atom + 1;
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

        int idx_cell = 0;
        int idx_atom = 0;

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
                    cal_fiJ(num_source, 0, cells, idx_cell, idx_atom);
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
        .def("cal_multipole", &cal_fmm::cal_multipole)
        .def("cal_M2M", &cal_fmm::cal_M2M)
        .def("cal_force", &cal_fmm::cal_force)
        .def_readwrite("natom", &cal_fmm::natom)
        .def_readwrite("n_crit", &cal_fmm::n_crit)
        .def_readwrite("theta", &cal_fmm::theta)
        .def_readwrite("charges", &cal_fmm::charges)
        .def_readwrite("t_crd", &cal_fmm::t_crd)
        .def_readwrite("gpair_table_fmm", &cal_fmm::gpair_table_fmm)
        .def_readwrite("gname_iatoms_pairs", &cal_fmm::gname_iatoms_pairs)
        ;
    
    py::class_<cal_fmm::Atomwise>(m, "Atomwise")
        .def(py::init<>())
        .def_readwrite("i", &cal_fmm::Atomwise::i)
        .def_readwrite("j", &cal_fmm::Atomwise::j)
        .def_readwrite("f", &cal_fmm::Atomwise::f)
        .def_readwrite("r", &cal_fmm::Atomwise::r)
        ;
    
    py::class_<cal_fmm::Cellwise>(m, "Cellwise")
        .def(py::init<>())
        .def_readwrite("atom_i", &cal_fmm::Cellwise::atom_i)
        .def_readwrite("atoms_J", &cal_fmm::Cellwise::atoms_J)
        .def_readwrite("f", &cal_fmm::Cellwise::f)
        .def_readwrite("r", &cal_fmm::Cellwise::r)
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