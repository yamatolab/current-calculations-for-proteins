#include<iostream>
#include<Eigen/Dense>
#include<cmath>
#include<pybind11/pybind11.h>
using namespace py = pybind11;
using namespace Eigen;

class cal_fmm{


    // public variables
    public:
        int n_crit;
        float theta;
        MatrixXd q;
        MatrixXd t_crd;
        MatrixXd t_pbc;
    

    void setup(int n_crit, float theta, MatrixXd q){
        n_crit = n_crit;
        theta = theta;
        q = q;
    };

    // read trajectory
    void initialize(const MatrixXd& crd, const MatrixXd & pbc){
        t_crd = crd;
        t_pbc = pbc;
        
    };

    struct Cell {
        int              nleaf;
        Eigen::VectorXi  leaf;
        int              nchild;
        Eigen::VectorXi  child;
        int              parent;
        Eigen::Vector3d  rc;
        double           r;
        Eigen::VectorXd  multipole;

        Cell(int n_crit = 10)
            : nleaf(0),                             // number of atoms(leaf) in the cell
            leaf(Eigen::VectorXi::Zero(n_crit)),    // index of atoms in the cell
            nchild(0),                              // number of child cells
            child(Eigen::VectorXi::Zero(8)),        // index of 8 child cells
            parent(0),                              // index of parent cell
            rc(Eigen::Vector3d::Zero()),            // center of the cell
            r(0.0),
            multipole(Eigen::VectorXd::Zero(10))    // 10 multipoles
        {}
    };

////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // max-min
    float calculate_rc_cand1(crd, particles){
        Vectorxd rc(3);
        Vectorxd particles = particles;
        for (int i = 0, i < 3, i++){
            float r_max;
            float r_min;
            float r_cand;

            for (int j = 0, j < particles.size(), j++){
                r_cand = t_crd(particles(j)-1, i);
                if (j == 0){
                    r_max = r_cand;
                    r_min = r_cand;
                }
                else{
                    if (r_cand > r_max){
                        r_max = r_cand;
                    }
                    if (r_cand < r_min){
                        r_min = r_cand;
                    }
                }
            }

            rc(i) = r_min + abs(r_max - r_min) * 0.5;
        }
        return rc;
    }
    // average
    float calculate_rc_cand2(crd, particles){
        Matrixmd rc(3);
        Vectorxd particles = particles;
        float inv_n = 1.0 / particles.size();

        for (int i = 0, i < 3, i++){
            Vectorxd t_coord = t_crd.col(i)
            float r_sum = 0.0; 

            for (int j = 0, j < particles.size(), j++){
                float r_sum = r_sum + t_coord(particles(j)-1);
            }

            rc(i) = r_sum * inv_n;
        }
        return rc;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////

    float cal_multipole(multipole, rc, particles){
        
        for (int i = 0, i < particles.size(), i++){
            float dx = rc(0) - t_crd(particles(i)-1, 0);
            float dy = rc(1) - t_crd(particles(i)-1, 1);
            float dz = rc(2) - t_crd(particles(i)-1, 2);
            float qj = 

            multipole(0) = multipole(0) + 1.0;
            multipole(1) = multipole(1) + dx;
            multipole(2) = multipole(2) + dy;
            multipole(3) = multipole(3) + dz;
            multipole(4) = multipole(4) + dx * dx;
            multipole(5) = multipole(5) + dy * dy;
            multipole(6) = multipole(6) + dz * dz;
            multipole(7) = multipole(7) + 2 * dx * dy;
            multipole(8) = multipole(8) + 2 * dy * dz;
            multipole(9) = multipole(9) + 2 * dz * dx;
            
        }
    }

    float cal_M2M(cell){
        Vectorxd p_potential(10);
        Vectorxd c_potential(10);
        Vectorxd c_rc(3);
        Vectorxd p_rc(3);
        

        // cellをよみこむ, cellの内容を読み込む
        dx = p_rc(0) - c_rc(0);
        dy = p_rc(1) - c_rc(1);
        dz = p_rc(2) - c_rc(2);

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

    float cal_fiJ(source, rc, potential, theta, cell){
        float crd_target = crd(target-1);
        float rx = crd_target(0) - rc(0);
        float ry = crd_target(1) - rc(1);
        float rz = crd_target(2) - rc(2);
        float r = sqrt(pow(rx, 2) + pow(ry, 2) + pow(rz, 2));

        if (cell.r > theta * r){
            bool do_fmm = False;
            exit();
        }
        else{

            float dx = crd_target(0) - rc(0);
            float dy = crd_target(1) - rc(1);
            float dz = crd_target(2) - rc(2);

            vector::Vectorxd bJx(10);
            vector::Vectorxd bJy(10);
            vector::Vectorxd bJz(10);
            float inv_r = 1.0 / r;
            float r2 = inv_r * inv_r;
            float r3 = r2 * inv_r;
            float r5 = r3 * r2;
            float r7 = r5 * r2;

            float dx2 = dx * dx;
            float dy2 = dy * dy;
            float dz2 = dz * dz;

            float dxdy = dx * dy;
            float dydz = dy * dz;
            float dzdx = dz * dx;

            float dxr5 = 3 * dx * r5;
            float dyr5 = 3 * dy * r5;
            float dzr5 = 3 * dz * r5;

            float dxdydz = 15 * dxdy * dz * r7;
            float dx2dy  = 15 * dx2 * dy * r7;
            float dy2dz  = 15 * dy2 * dz * r7;
            float dz2dx  = 15 * dz2 * dx * r7;
            float dy2dx  = 15 * dy2 * dx * r7;
            float dz2dy  = 15 * dz2 * dy * r7;
            float dx2dz  = 15 * dx2 * dz * r7;


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
            float fx = potential * bJx;
            float fy = potential * bJy;
            float fz = potential * bJz;

            return vector::Vectorxd(fx, fy, fz), vector::Vectorxd(rx, ry, rz);
        }
    }
}




py(lib_fmm, m){
    
}