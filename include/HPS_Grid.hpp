#ifndef HPS_GRID_HPP_
#define HPS_GRID_HPP_

#include "HPS_Vector.hpp"

// @@TODO: Add documentation to all

namespace hps {

enum DIMS {
    X,
    Y,
    Z
};

template<class T, int D>
class Grid {

public:

    Grid() {}
    virtual T point(int dim, int index) = 0;

};

template<class T, int D>
class EdgeGrid : public Grid<T, D> {

public:

    int N_pts[D];
    T lower_limit[D];
    T upper_limit[D];
    hps::Vector<T> data_[D];
    T spacing[D];

    EdgeGrid(int N_points[], T lower[], T upper[]) {
        for (int d = 0; d < D; d++) {
            this->N_pts[d] = N_points[d];
            this->lower_limit[d] = lower[d];
            this->upper_limit[d] = upper[d];
            this->data_[d] = hps::Vector<T>(N_points[d]);
            this->spacing[d] = (upper[d] - lower[d]) / (N_points[d] - 1);
        }

        for (int d = 0; d < D; d++) {
            for (int i = 0; i < this->N_pts[d]; i++) {
                this->data_[d][i] = lower[d] + i*this->spacing[d];
            }
        }
    }

    T point(int dim, int index) {
        return this->data_[dim][index];
    }

};

template<class T, int D>
class CellGrid : public Grid<T, D> {

public:

    int N_pts[D];
    T lower_limit[D];
    T upper_limit[D];
    hps::Vector<T> data_[D];
    T spacing[D];
    
    CellGrid(int N_points, T lower, T upper) {
        this->N_pts[0] = N_points;
        this->lower_limit[0] = lower;
        this->upper_limit[0] = upper;
        this->data_[0] = hps::Vector<T>(N_points);
        this->spacing[0] = (upper - lower) / (N_points);
        
        for (int i = 0; i < N_points; i++) {
            this->data_[0][i] = (lower + this->spacing[0]/2) + i*this->spacing[0];
        }
    }

    CellGrid(int N_points[], T lower[], T upper[]) {
        for (int d = 0; d < D; d++) {
            this->N_pts[d] = N_points[d];
            this->lower_limit[d] = lower[d];
            this->upper_limit[d] = upper[d];
            this->data_[d] = hps::Vector<T>(N_points[d]);
            this->spacing[d] = (upper[d] - lower[d]) / (N_points[d]);
        }

        for (int d = 0; d < D; d++) {
            for (int i = 0; i < this->N_pts[d]; i++) {
                this->data_[d][i] = (lower[d] + this->spacing[d]/2) + i*this->spacing[d];
            }
        }
    }

    T point(int dim, int index) {
        return this->data_[dim][index];
    }

};

}

#endif // HPS_GRID_HPP_

