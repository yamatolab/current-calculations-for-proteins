#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
namespace py = pybind11;

class get_table_fmm {

public:
    // public variables
    std::vector<std::pair<int, int>> bonded_pairs;
    std::vector<std::pair<std::string, std::vector<std::string>>> gpair_table;
    std::vector<std::pair<std::string, std::vector<int>>> gname_iatoms_pairs;

    std::vector<std::pair<int, int>> extracted_bonded_pairs;

    get_table_fmm() :
        bonded_pairs(std::vector<std::pair<int, int>>()),
        gpair_table(std::vector<std::pair<std::string, std::vector<std::string>>>()),
        gname_iatoms_pairs(std::vector<std::pair<std::string, std::vector<int>>>()),
        extracted_bonded_pairs(std::vector<std::pair<int, int>>())
    {};


    void setup(const std::vector<std::pair<int, int>>& input_bonded_pairs,
        const std::vector<std::pair<std::string, std::vector<int>>>& input_gname_iatoms_pairs,
        const std::vector<std::pair<std::string, std::vector<std::string>>>& input_gpair_table) {
        
        bonded_pairs = input_bonded_pairs;
        gname_iatoms_pairs = input_gname_iatoms_pairs;
        gpair_table = input_gpair_table;
    };

    // extract bonded pairs based on group pairs
    std::vector<std::pair<int, int>> extract_bonded_pairs(){

        int size_gpair_table = gpair_table.size();
        int size_bonded_pairs = bonded_pairs.size();

        for (int i = 0; i < size_gpair_table; i++) {
            std::string gname_i = gpair_table[i].first;
            std::vector<std::string> gname_js = gpair_table[i].second;
            int size_gname_js = gname_js.size();

            for (int j = 0; j < size_gname_js; j++) {
                std::string gname_j = gname_js[j];
                std::vector<int> atoms_i = get_atoms(gname_i, gname_iatoms_pairs);
                std::vector<int> atoms_j = get_atoms(gname_j, gname_iatoms_pairs);

                for (int k = 0; k < size_bonded_pairs; k++) {
                    int atom1 = bonded_pairs[k].first;
                    int atom2 = bonded_pairs[k].second;

                    // Check if atom1 is in group i and atom2 is in group j
                    if (std::find(atoms_i.begin(), atoms_i.end(), atom1) != atoms_i.end() &&
                        std::find(atoms_j.begin(), atoms_j.end(), atom2) != atoms_j.end()) {
                        extracted_bonded_pairs.push_back({ atom1, atom2 });
                    }
                }
            }
        }
        return extracted_bonded_pairs;
    };

    std::vector<int> get_atoms(const std::string& gname, const std::vector<std::pair<std::string, std::vector<int>>>& gname_iatoms_pairs) {

        std::vector<int> atoms;
        for (const auto& pair : gname_iatoms_pairs) {
            if (pair.first == gname) {
                atoms = pair.second;
                break;
            }
        }
        return atoms;
    };

    // std::vector<std::pair<int, std::vector<int>>> get_gnum_iatoms_pairs() {

    //     int size_gname_iatoms_pairs = gname_iatoms_pairs.size();

    //     for (int i = 0; i < size_gname_iatoms_pairs; i++) {

    //         std::vector<int> iatoms = gname_iatoms_pairs[i].second;
    //         int gnum = i + 1;;
    //         gnum_iatoms_pairs.push_back({ num_iatoms, iatoms });
    //     }
    //     return gnum_iatoms_pairs
    // };
    
};

PYBIND11_MODULE(lib_table_fmm, m) {
    py::class_<get_table_fmm>(m, "get_table_fmm")
        .def(py::init<>())
        .def("setup", &get_table_fmm::setup)
        .def("extract_bonded_pairs", &get_table_fmm::extract_bonded_pairs)
        .def_readwrite("bonded_pairs", &get_table_fmm::bonded_pairs)
        .def_readwrite("gpair_table", &get_table_fmm::gpair_table)
        .def_readwrite("gname_iatoms_pairs", &get_table_fmm::gname_iatoms_pairs)
        .def_readwrite("extracted_bonded_pairs", &get_table_fmm::extracted_bonded_pairs)
        ;
}
