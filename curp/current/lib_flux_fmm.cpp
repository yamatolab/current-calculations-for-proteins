#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;
using namespace Eigen;



class cal_flux_fmm{

    public:
        int natom;
        int ngrp;
        std::vector<int> target_atoms;
        std::vector<int> iatoms_group;
        MatrixXd t_vel;
        MatrixXd hflux_ij;
        MatrixXd eflux_ij;
        bool flag_heat = false;
        bool flag_energy = false;

        struct Atomwise{
            int atom_i;
            int atom_j;
            Vector3d f;
            Vector3d r;

        };

        struct Cellwise{
            int atom_i;
            VectorXi atoms_J;
            Vector3d f;
            Vector3d r;

        };

        std::vector<Atomwise> atomwise;
        std::vector<Cellwise> cellwise;


    void initialize(std::vector<int>& target_atoms_input, std::vector<int>& iatoms_group_input){
        target_atoms = target_atoms_input;
        iatoms_group = iatoms_group_input;
        natom = target_atoms.size();
        ngrp = iatoms_group.size();

    };

    void set_flux(const std::string& flux_type){

        if (flux_type == "heat"){
            flag_heat = true;
            hflux_ij = MatrixXd::Zero(ngrp, ngrp);
        }
        else if (flux_type == "energy"){
            flag_energy = true;
            eflux_ij = MatrixXd::Zero(ngrp, ngrp);
        }
        else{
            std::cout << "Invalid flux type" << std::endl;
        }
    };

    void init_cal(const MatrixXd& vel){

        t_vel = vel;
        if (flag_heat == true){
            hflux_ij = MatrixXd::Zero(ngrp, ngrp);
        }
        else if (flag_energy == true){
            eflux_ij = MatrixXd::Zero(ngrp, ngrp);
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //heat flux//

    void cal_hflux_atomwise(const std::vector<Atomwise>& atomwise){

        int size_atomwise = atomwise.size();
        for (int i = 0; i < size_atomwise; i++){
            int atom_i = atomwise[i].atom_i;
            int atom_j = atomwise[i].atom_j;
            int idx_i = atom_i - 1;
            int idx_j = atom_j - 1;

            Vector3d fij = atomwise[i].f;
            Vector3d r = atomwise[i].r;

            Vector3d vi = t_vel.col(idx_i);
            Vector3d vj = t_vel.col(idx_j);
            Vector3d vij = vi + vj;

            Vector3d h_ij = r * (fij.dot(vij)) * 0.5;

            int igrp = iatoms_group[idx_i];
            int jgrp = iatoms_group[idx_j];

            if (igrp == 0 or jgrp == 0){
                continue;
            }
            
            hflux_ij(igrp, jgrp) += h_ij;        
        }
    };

    void cal_hflux_cellwise(std::vector<Cellwise> cellwise){

        int size_cellwise = cellwise.size();
        for (int i = 0; i < size_cellwise; i++){
            int atom_i = cellwise[i].atom_i;
            int idx_i = atom_i - 1;

            Vector3d fij = cellwise[i].f;
            Vector3d r = cellwise[i].r;

            Vector3d vi = t_vel.col(idx_i);

            Vector3d h_ij = r * (fij.dot(vi)) * 0.5;

            int atom_J = cellwise[i].atoms_J.coeff(0);
            int igrp = iatoms_group[idx_i];
            int Jgrp = iatoms_group[atom_J - 1];

            if (igrp == 0 or Jgrp == 0){
                continue;
            }

            hflux_ij(igrp, Jgrp) += h_ij;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //energy flux//

    void cal_eflux_atomwise(const std::vector<Atomwise>& atomwise){

        int size_atomwise = atomwise.size();
        for (int i = 0; i < atomwise.size(); i++){
            int atom_i = atomwise[i].atom_i;
            int atom_j = atomwise[i].atom_j;
            int idx_i = atom_i - 1;
            int idx_j = atom_j - 1;

            Vector3d fij = atomwise[i].f;
            Vector3d r = atomwise[i].r;

            Vector3d vi = t_vel.col(idx_i);
            Vector3d vj = t_vel.col(idx_j);
            Vector3d vij = vi + vj;

            double e_ij = fij.dot(vij);

            int igrp = iatoms_group[idx_i];
            int jgrp = iatoms_group[idx_j];

            if (igrp == 0 or jgrp == 0){
                continue;
            }

            eflux_ij(igrp, jgrp) += e_ij;

        }
    };

    void cal_eflux_cellwise(std::vector<Cellwise> cellwise){

        int size_cellwise = cellwise.size();
        for (int i = 0; i < cellwise.size(); i++){
            int atom_i = cellwise[i].atom_i;
            int idx_i = atom_i - 1;

            Vector3d fij = cellwise[i].f;
            Vector3d r = cellwise[i].r;

            Vector3d vi = t_vel.col(idx_i);

            double e_ij = fij.dot(vi);

            int atom_J = cellwise[i].atoms_J(0);
            int igrp = iatoms_group[idx_i];
            int Jgrp = iatoms_group[atom_J - 1];
            
            if (igrp == 0 or jgrp == 0){
                continue;
            }

            eflux_ij(igrp, Jgrp) += e_ij;

        }
    };
}

PYBIND11_MODULE(lib_flux_fmm, m){
    py::class_<cal_flux_fmm>(m, "cal_flux_fmm")
        .def(py::init<>())
        .def("initialize", &cal_flux_fmm::initialize)
        .def("set_flux", &cal_flux_fmm::set_flux)
        .def("init_cal", &cal_flux_fmm::init_cal)
        .def("cal_hflux_atomwise", &cal_flux_fmm::cal_hflux_atomwise)
        .def("cal_hflux_cellwise", &cal_flux_fmm::cal_hflux_cellwise)
        .def("cal_eflux_atomwise", &cal_flux_fmm::cal_eflux_atomwise)
        .def("cal_eflux_cellwise", &cal_flux_fmm::cal_eflux_cellwise)
        .def_readwrite("natom", &cal_flux_fmm::natom)
        .def_readwrite("ngrp", &cal_flux_fmm::ngrp)
        .def_readwrite("target_atoms", &cal_flux_fmm::target_atoms)
        .def_readwrite("iatoms_group", &cal_flux_fmm::iatoms_group)
        .def_readwrite("t_vel", &cal_flux_fmm::t_vel)
        .def_readwrite("atomwise", &cal_flux_fmm::atomwise)
        .def_readwrite("cellwise", &cal_flux_fmm::cellwise)
        .def_readwrite("flag_heat", &cal_flux_fmm::flag_heat)
        .def_readwrite("flag_energy", &cal_flux_fmm::flag_energy)
        .def_readwrite("hflux_ij", &cal_flux_fmm::hflux_ij)
        .def_readwrite("eflux_ij", &cal_flux_fmm::eflux_ij);
}