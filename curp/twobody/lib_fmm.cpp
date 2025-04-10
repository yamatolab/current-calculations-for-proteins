#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;
using namespace Eigen;

class cal_fmm{


    // public variables
    public:
        int natom;
        int n_crit;
        double theta;
        VectorXd charges;
        MatrixXd t_crd;

        struct Gpair_table_fmm{
            std::string group_i;
            std::vector<std::string> group_js;
            
            Gpair_table_fmm():
                group_i(""),
                group_js(std::vector<std::string>())
            {}
        };

        std::vector<Gpair_table_fmm> gpair_table_fmm;

        struct Gname_iatoms_pairs{
            std::string group;
            VectorXi iatoms;

            Gname_iatoms_pairs():
                group(""),
                iatoms(VectorXi::Zero(0))
            {}
        };

        std::vector<Gname_iatoms_pairs> gname_iatoms_pairs;

        struct Atomwise{
            int i;
            int j;
            Vector3d f;
            Vector3d r;
    
            Atomwise():
                i(0),
                j(0),
                f(Vector3d::Zero()),
                r(Vector3d::Zero())
            {}
        };

        std::vector<Atomwise> atomwise;
    
        struct Cellwise{
            int atom_i;         //atom number of the source(i)
            VectorXi atoms_J;   //atom numbers of the target cell J
            Vector3d f;         //force between atom i and Cell J
            Vector3d r;         //distance between atom i and center of Cell J
    
            Cellwise():
                atom_i(0),
                atoms_J(VectorXi::Zero(0)),
                f(Vector3d::Zero()),
                r(Vector3d::Zero())
            {}
        };

        std::vector<Cellwise> cellwise;

        int idx_cell;
        int idx_atom;

    void setup(const int input_natom, const int& input_n_crit, const double& input_theta, \
        const VectorXd& input_charges, std::vector<Gname_iatoms_pairs>& input_gname_iatoms_pairs){
        natom  = input_natom;
        n_crit = input_n_crit;
        theta  = input_theta;
        charges = input_charges;
        gname_iatoms_pairs = input_gname_iatoms_pairs;
    };

    // read trajectory
    void initialize(const MatrixXd& crd){
        t_crd = crd;
        atomwise = std::vector<Atomwise>();
        cellwise = std::vector<Cellwise>();

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // calculate the center and radius of the cell 
    void calculate_rc(std::vector<int> atoms){

        Vector3d rc = Vector3d::Zero();
        Vector3d r = Vector3d::Zero();
        int size_atoms = atoms.size();

        MatrixXd crd_atom = VectorXd::Zero(size_atoms,3);

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
    };
        
    

////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void cal_multipole(VectorXd& multipole, Vector3d& rc, VectorXi& atoms){
        
        for (int i = 0; i < std::size(atoms); i++){

            int index_atom = atoms(i) - 1;

            Vector3d crd_atom = t_crd.row(index_atom);
            double dx = rc(0) - crd_atom(0);
            double dy = rc(1) - crd_atom(1);
            double dz = rc(2) - crd_atom(2);
            double qj = charges(index_atom);

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
    };

    void cal_M2M(std::vector<Cell> cells){

        int cells_size = cells.size();
        for (int i = 0; i < cells_size; i++){
            int i_inv = cells.size() - 1 - i;

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
                        double charge = charges(idx_source);
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

                double charge_i = charges(idx_source);
                double charge_j = charges(idx_target);
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

    std::vector<int> get_atoms(const std::string& source, const std::vector<Gname_iatoms_pairs>& pairs) {
        for (const auto& pair : pairs) {
            if (pair.group == source) {
                return std::vector<int>(pair.iatoms.data(), pair.iatoms.data() + pair.iatoms.size());
            }
        }
        return std::vector<int>();
    };

    std::vector<Cell> get_cells(const std::string& source, const std::vector<All_cells>& all_cells) {
        for (const auto& all_cell : all_cells) {
            if (all_cell.group == source) {
                return all_cell.cells;
            }
        }
        return std::vector<Cell>();
    };

    void cal_force(const std::vector<All_cells>& all_cells){

        int idx_cell = 0;
        int idx_atom = 0;

        int size_table = gpair_table_fmm.size();
        for (int i = 0; i < size_table; i++){

            std::string source = gpair_table_fmm[i].group_i;
            std::vector<std::string> targets = gpair_table_fmm[i].group_js;
        
            std::vector<int> source_atoms = get_atoms(source, gname_iatoms_pairs);
            int source_size = source_atoms.size();

            for (int i = 0; i < source_size; i++){

                int num_source = source_atoms[i];
                int size_targets = targets.size();

                for (int j = 0; j < size_targets; j++){      

                    std::string target = targets[j];
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
        .def_readwrite("idx_cell", &cal_fmm::idx_cell)
        .def_readwrite("idx_atom", &cal_fmm::idx_atom)
        ;
    
}