#include "HPS_Patch.hpp"

namespace hps {

Patch::Patch() {}
Patch::Patch(hps::CellGrid<double, 2> patch_grid, int ID, int level, bool is_leaf) : grid(patch_grid), ID(ID), level(level), is_leaf(is_leaf) {}

std::pair<Patch, Patch> Patch::split(int first_ID, int second_ID) {
    if (this->is_leaf) {
        if (this->level % 2 == 0) {
            // Even level, split square patch into two rectangles
            // Create two new grids for the new patches from the original grid
            //    Lower grid
            double lower_grid_lower_limits[2] = {this->grid.lower_limit[X], this->grid.lower_limit[Y]};
            double lower_grid_upper_limits[2] = {this->grid.upper_limit[X], this->grid.upper_limit[Y]/2};
            CellGrid<double, 2> lower_grid(this->grid.N_pts, lower_grid_lower_limits, lower_grid_upper_limits);

            //    Upper grid
            double upper_grid_lower_limits[2] = {this->grid.lower_limit[X], this->grid.upper_limit[Y]/2};
            double upper_grid_upper_limits[2] = {this->grid.upper_limit[X], this->grid.upper_limit[Y]};
            CellGrid<double, 2> upper_grid(this->grid.N_pts, upper_grid_lower_limits, upper_grid_upper_limits);

            // Create new patches metainfo
            this->is_leaf = false;
            Patch lower_patch(lower_grid, first_ID, level+1, true);
            Patch upper_patch(upper_grid, second_ID, level+1, true);
            return std::pair<Patch, Patch>(lower_patch, upper_patch);
        }
        else {
            // Odd level, split rectangle patch into two squares
            // Create two new grids for the new patches from the original grid
            //    Left grid
            double left_grid_lower_limits[2] = {this->grid.lower_limit[X], this->grid.lower_limit[Y]};
            double left_grid_upper_limits[2] = {this->grid.upper_limit[X]/2, this->grid.upper_limit[Y]};
            CellGrid<double, 2> left_grid(this->grid.N_pts, left_grid_lower_limits, left_grid_upper_limits);

            //    Right grid
            double right_grid_lower_limits[2] = {this->grid.upper_limit[X]/2, this->grid.lower_limit[Y]};
            double right_grid_upper_limits[2] = {this->grid.upper_limit[X], this->grid.upper_limit[Y]};
            CellGrid<double, 2> right_grid(this->grid.N_pts, right_grid_lower_limits, right_grid_upper_limits);

            // Create new patches metainfo
            this->is_leaf = false;
            Patch left_patch(left_grid, first_ID, level+1, true);
            Patch right_patch(right_grid, second_ID, level+1, true);
            return std::pair<Patch, Patch>(left_patch, right_patch);
        }
    }
    else {
        throw std::invalid_argument("[Patch::split] Patch being split is not a leaf");
    }
}

Matrix<double> Patch::buildDirichletToNeumannMatrix() {
    if (!this->is_leaf) {
        throw std::invalid_argument("[Patch::buildDirichletToNeumannMatrix] Patch to build DtN matric for is not a leaf");
    }

    return buildDirichletToNeumann(this->grid);
}



}
