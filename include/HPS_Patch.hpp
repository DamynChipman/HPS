#ifndef HPS_PATCH_HPP_
#define HPS_PATCH_HPP_

#include "HPS_Base.hpp"
#include "HPS_PatchSolver.hpp"
#include <utility>

namespace hps {

class Patch {

public:

    hps::CellGrid<double, 2> grid;
    int ID;
    int level;
    bool is_leaf;

    Patch();
    Patch(hps::CellGrid<double, 2> patch_grid, int ID, int level, bool is_leaf);
    std::pair<Patch, Patch> split(int first_ID, int second_ID);
    Matrix<double> buildDirichletToNeumannMatrix(double lambda);

};

}

#endif // HPS_PATCH_HPP_
