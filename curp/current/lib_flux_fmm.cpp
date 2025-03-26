#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
using namespace py = pybind11;
using namespace Eigen;



class cal_flux_fmm{

    public:
        int natom;
        VectorXd target_atoms;
        VectorXd iatoms_group;
        MatrixXd t_vel;
        MatrixXd hflux_ij;
        MatrixXd eflux_ij;


        struct Atomwise{
            int atom_i;
            int atom_j;
            VectorXd f;
            VectorXd r;
            
        }

        struct Cellwise{
            std::string group_i;
            std::string group_j;
            int atom_i;
            VectorXd f;
            VectorXd r;

        };

        std::vector<Atomwise> atomwise;
        std::vector<Cellwise> cellwise;


    void initialize(std::vector<int>& target_atoms_input, std::vector<int>& iatoms_group_input){
        target_atoms = target_atoms_input;
        iatoms_group = iatoms_group_input;
        natom = target_atoms.size();

    }

    void set_flux(std::string flux_type){
        if (flux_type == "heat"){
            hflux_ij = MatrixXd::Zero(natom, natom);
        }
        else if (flux_type == "energy"){
            eflux_ij = MatrixXd::Zero(natom, natom);
        }
        else{
            std::cout << "Invalid flux type" << std::endl;
        }
    }

    void init_cal(MatrixXd& vel){
        t_vel = vel;
    }

    void cal_hflux_atomwise(std::vector<Atomwise>& atomwise){

        for (int i = 0; i < atomwise.size(); i++){
            int atom_i = atomwise[i].atom_i;
            int atom_j = atomwise[i].atom_j;
            Vector3d fij = atomwise[i].f;
            Vector3d r = atomwise[i].r;

            Vector3d vi = t_vel.col(atom_i - 1);
            Vector3d vj = t_vel.col(atom_j - 1);
            Vector3d vij = vi + vj;

            Vector3d hflux_ij = r * (fij.dot(vij));
            
            igrp = iatoms_group[atom_i - 1];
            jgrp = iatoms_group[atom_j - 1];
            hflux[igrp, jgrp] += hflux_ij;

        }
    }
        
    void cal_hflux_cellwise(std::vector<Cellwise> cellwise){
        for (int i = 0; i < cellwise.size(); i++){
            string group_i = cellwise[i].group_i;
            string group_j = cellwise[i].group_j;
            int atom_i = cellwise[i].atom_i;
            Vector3d fij = cellwise[i].f;
            Vector3d r = cellwise[i].r;

            Vector3d vi = t_vel.col(atom_i - 1);

            Vector3d hflux_ij = r * (fij.dot(vij));
            
            hflux[group_i, group_j] += hflux_ij;
        }
    }
}