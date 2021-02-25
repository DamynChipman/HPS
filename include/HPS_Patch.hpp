#ifndef HPS_PATCH_HPP_
#define HPS_PATCH_HPP_

#include "HPS_Base.hpp"
#include "HPS_PatchSolver.hpp"
#include <utility>

namespace hps {

class Patch {

public:

    // Metadata
    int ID;
    int level;
    bool is_leaf;
    int coords[2];

    // Patch domain
    hps::CellGrid<double, 2> grid;
    // int Nx = grid.N_pts[X];
    // int Ny = grid.N_pts[Y];

    // Local poisson data
    Matrix<double> f{}; // = Matrix<double>(Nx, Ny);

    // Local solution
    Matrix<double> u{};

    // Patch boundary data
    Vector<double> g{};

    // HPS data
    // Matrix<double> T = Matrix<double>(2*Nx + 2*Ny, 2*Nx + 2*Ny);
    Matrix<double> T{}; // = Matrix<double>(1, 1);
    Matrix<double> S{}; // = Matrix<double>(1, 1);
    Vector<double> fhat{}; // = Vector<double>(2*Nx + 2*Ny);

    Patch();
    Patch(hps::CellGrid<double, 2> parent_grid);
    Patch(hps::CellGrid<double, 2> patch_grid, int ID, int level, bool is_leaf);
    std::pair<Patch, Patch> split(int first_ID, int second_ID);
    Matrix<double> buildDirichletToNeumannMatrix(double lambda);
    void buildT(double lambda);

};

}

#endif // HPS_PATCH_HPP_
