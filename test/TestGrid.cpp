#include "gtest/gtest.h"
#include "HPS_Grid.hpp"

/*///////////////////////////////////////////////
 *
 *  Unit Tests -- CellGrid Class
 *
 *///////////////////////////////////////////////

TEST(CellGrid, init_1D) {
    
    int N = 4;
    double A = 0;
    double B = 1;
    hps::CellGrid<double, 1> grid(N, A, B);
    
    ASSERT_EQ(grid.N_pts[hps::X], N);
    ASSERT_EQ(grid.lower_limit[hps::X], A);
    ASSERT_EQ(grid.upper_limit[hps::X], B);
    ASSERT_EQ(grid.spacing[hps::X], 0.25);
    
}

TEST(CellGrid, init_2D) {
    
    int Nx = 3;
    double Ax = -1;
    double Bx = 1;

    int Ny = 5;
    double Ay = 0;
    double By = 2;
    
    int N[2] = {Nx, Ny};
    double A[2] = {Ax, Ay};
    double B[2] = {Bx, By};
    hps::CellGrid<double, 2> grid(N, A, B);
    
    ASSERT_EQ(grid.N_pts[hps::X], Nx);
    ASSERT_EQ(grid.lower_limit[hps::X], Ax);
    ASSERT_EQ(grid.upper_limit[hps::X], Bx);
    ASSERT_EQ(grid.spacing[hps::X], (double) 2.0/3.0);
    
    ASSERT_EQ(grid.N_pts[hps::Y], Ny);
    ASSERT_EQ(grid.lower_limit[hps::Y], Ay);
    ASSERT_EQ(grid.upper_limit[hps::Y], By);
    ASSERT_EQ(grid.spacing[hps::Y], (double) 2.0/5.0);
    
}

TEST(CellGrid, points) {
    
    int N = 4;
    double A = 0;
    double B = 1;
    hps::CellGrid<double, 1> grid(N, A, B);
    
    double points_test[4] = {0.125, 0.375, 0.625, 0.875};
    
    for (int i = 0; i < N; i++) {
        ASSERT_FLOAT_EQ(grid.point(hps::X, i), points_test[i]);
    }
    
}