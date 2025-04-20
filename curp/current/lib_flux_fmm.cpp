#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
namespace py = pybind11;
using namespace Eigen;



class cal_flux_fmm{

public:
    int natom;
    int ngrp;
    std::vector<int> target_atoms;
    VectorXi iatoms_group;
    MatrixXd t_vel;
    std::vector<std::vector<Vector3d>> hflux_ij;
    MatrixXd eflux_ij;
    bool flag_heat = false;
    bool flag_energy = false;

    // create instance
    cal_flux_fmm():
        natom(0),
        ngrp(0),
        target_atoms(std::vector<int>()),
        iatoms_group(VectorXi::Zero(0)),
        t_vel(MatrixXd::Zero(natom, 3)),
        hflux_ij(std::vector<std::vector<Vector3d>>()),
        eflux_ij(MatrixXd::Zero(0, 0)),
        flag_heat(false),
        flag_energy(false)
    {};


    void initialize(std::vector<int>& target_atoms_input, VectorXi& iatoms_group_input){
        target_atoms = target_atoms_input;
        iatoms_group = iatoms_group_input;
        natom = target_atoms.size();
        ngrp = iatoms_group.size();

    };

    void set_flux(const std::string& flux_type){

        if (flux_type == "heat"){
            flag_heat = true;
            hflux_ij.resize(ngrp, std::vector<Vector3d>(ngrp, Vector3d::Zero()));
        }
        else if (flux_type == "energy"){
            flag_energy = true;
            eflux_ij = MatrixXd::Zero(ngrp, ngrp);
        }
        else{
            std::cerr << "Invalid flux type" << std::endl;
            std::exit(1);
        }
    };

    void init_cal(const MatrixXd& vel){

        t_vel = vel;
        if (flag_heat == true){
            hflux_ij.resize(natom, std::vector<Vector3d>(natom, Vector3d::Zero()));
        }
        else if (flag_energy == true){
            eflux_ij = MatrixXd::Zero(natom, natom);
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //heat flux//

    void cal_hflux_atomwise(const std::vector<int>& atomwise_i, const std::vector<int>& atomwise_j, \
                            const std::vector<Vector3d>& atomwise_f, const std::vector<Vector3d>& atomwise_r){

        int size_atomwise = atomwise_i.size();
        for (int i = 0; i < size_atomwise; i++){
            int atom_i = atomwise_i[i];
            int atom_j = atomwise_j[i];
            int idx_i = atom_i - 1;
            int idx_j = atom_j - 1;

            Vector3d fij = atomwise_f[i];
            Vector3d r = atomwise_r[i];

            Vector3d vi = t_vel.row(idx_i);
            Vector3d vj = t_vel.row(idx_j);
            Vector3d vij = vi + vj;

            Vector3d h_ij = r * (fij.dot(vij)) * 0.5;

            int igrp = iatoms_group(idx_i);
            int jgrp = iatoms_group(idx_j);

            if (igrp == 0 or jgrp == 0){
                continue;
            }
            
            if (igrp != jgrp){
                hflux_ij[igrp-1][jgrp-1] += h_ij;
                hflux_ij[jgrp-1][igrp-1] += h_ij;
            }
            else{
                hflux_ij[igrp-1][jgrp-1] += 0.5 * h_ij;
            }
        }
    };

    void cal_hflux_cellwise(const std::vector<int>& cellwise_i, const std::vector<std::vector<int>>& cellwise_J, \
                            const std::vector<Vector3d>& cellwise_f, const std::vector<Vector3d>& cellwise_r){

        int size_cellwise = cellwise_i.size();
        for (int i = 0; i < size_cellwise; i++){
            int atom_i = cellwise_i[i];
            int idx_i = atom_i - 1;

            Vector3d fij = cellwise_f[i];
            Vector3d r = cellwise_r[i];

            Vector3d vi = t_vel.row(idx_i);

            Vector3d h_ij = r * (fij.dot(vi)) * 0.5;

            int atom_J = cellwise_J[i][0];
            int igrp = iatoms_group(idx_i);
            int Jgrp = iatoms_group(atom_J - 1);

            if (igrp == 0 or Jgrp == 0){
                continue;
            }

            if (igrp != Jgrp){
                hflux_ij[igrp-1][Jgrp-1] += h_ij;
                hflux_ij[Jgrp-1][igrp-1] += h_ij;
            }
            else{
                hflux_ij[igrp-1][Jgrp-1] += 0.5 * h_ij;
            }
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //energy flux//

    void cal_eflux_atomwise(std::vector<int>& atomwise_i, std::vector<int>& atomwise_j, \
                            std::vector<Vector3d>& atomwise_f){

        int size_atomwise = atomwise_i.size();
        for (int i = 0; i < size_atomwise; i++){
            int atom_i = atomwise_i[i];
            int atom_j = atomwise_j[i];
            int idx_i = atom_i - 1;
            int idx_j = atom_j - 1;

            Vector3d fij = atomwise_f[i];

            Vector3d vi = t_vel.row(idx_i);
            Vector3d vj = t_vel.row(idx_j);
            Vector3d vij = vi + vj;

            double e_ij = fij.dot(vij);

            int igrp = iatoms_group(idx_i);
            int jgrp = iatoms_group(idx_j);

            if (igrp == 0 or jgrp == 0){
                continue;
            }

            if (igrp != jgrp){
                eflux_ij(igrp-1, jgrp-1) += e_ij;
                eflux_ij(jgrp-1, igrp-1) += -e_ij;
            }
            else{
                eflux_ij(igrp-1, jgrp-1) += 0.5 * e_ij;
            }

        }
    };

    void cal_eflux_cellwise(std::vector<int>& cellwise_i, std::vector<std::vector<int>>& cellwise_J, \
                            std::vector<Vector3d>& cellwise_f){

        int size_cellwise = cellwise_i.size();
        for (int i = 0; i < size_cellwise; i++){
            int atom_i = cellwise_i[i];
            int idx_i = atom_i - 1;

            Vector3d fij = cellwise_f[i];

            Vector3d vi = t_vel.row(idx_i);

            double e_ij = fij.dot(vi);

            int atom_J = cellwise_J[i][0];
            int igrp = iatoms_group(idx_i);
            int Jgrp = iatoms_group(atom_J - 1);
            
            if (igrp == 0 or Jgrp == 0){
                continue;
            }

            if (igrp != Jgrp){
                eflux_ij(igrp-1, Jgrp-1) += e_ij;
                eflux_ij(Jgrp-1, igrp-1) += -e_ij;
            }
            else{
                eflux_ij(igrp-1, Jgrp-1) += 0.5 * e_ij;
            }
        }
    };
};

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
        .def_readwrite("flag_heat", &cal_flux_fmm::flag_heat)
        .def_readwrite("flag_energy", &cal_flux_fmm::flag_energy)
        .def_readwrite("hflux_ij", &cal_flux_fmm::hflux_ij)
        .def_readwrite("eflux_ij", &cal_flux_fmm::eflux_ij);

};
