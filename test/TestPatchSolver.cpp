#include "gtest/gtest.h"
#include "HPS_Base.hpp"
#include "HPS_PatchSolver.hpp"

/*///////////////////////////////////////////////
 *
 *  Unit Tests -- PatchSolver Functions
 *
 *///////////////////////////////////////////////

TEST(PatchSolver, mapSolution) {


    // Setup vectors with iteration info
    int N_runs = 10;
    hps::Vector<int> N_cells(N_runs);
    hps::Vector<double> u_errors(N_runs);
    hps::Vector<double> h_errors(N_runs);
    for (int i = 0; i < N_runs; i++) {
        N_cells[i] = pow(2, i+2);
    }

    // Begin convergence analysis loop
    for (int i = 0; i < N_runs; i++) {

        // Setup problem inputs
        int N_cell = N_cells[i];
        int size = N_cell * N_cell;
        double Ax = -2.0;
        double Bx = 1.0;
        double Ay = -3.0;
        double By = 5.0;
        int Nx = N_cell;
        int Ny = N_cell;
        int prob_ID = 0;

        // Create finite difference grid
        int N_pts[2] = {Nx, Ny};
        double lower_bounds[2] = {Ax, Ay};
        double upper_bounds[2] = {Bx, By};
        hps::CellGrid<double, 2> grid(N_pts, lower_bounds, upper_bounds);

        // Create Poisson problem
        hps::PoissonProblem poisson(prob_ID, Ax, Bx, Ay, By);

        // Build grid function vectors
        hps::Vector<double> g0(N_cell);              // East Dirichlet data
        hps::Vector<double> g1(N_cell);              // West Dirichlet data
        hps::Vector<double> g2(N_cell);              // South Dirichlet data
        hps::Vector<double> g3(N_cell);              // North Dirichlet data
        hps::Vector<double> g(4*N_cell);             // Stacked Dirichlet data
        hps::Vector<double> h0(N_cell);              // East Neumann data
        hps::Vector<double> h1(N_cell);              // West Neumann data
        hps::Vector<double> h2(N_cell);              // South Neumann data
        hps::Vector<double> h3(N_cell);              // North Neumann data
        hps::Vector<double> h_exact(4*N_cell);             // Stacked Neumann data
        hps::Matrix<double> f(N_cell, N_cell);       // RHS matrix
        hps::Matrix<double> u_exact(N_cell, N_cell); // Exact solution matrix

        // Fill in X data
        for (int i = 0; i < N_cell; i++) {
            double x = grid.point(hps::X, i);
            g2[i] = poisson.gSouth(x);
            g3[i] = poisson.gNorth(x);
            h2[i] = poisson.hSouth(x);
            h3[i] = poisson.hNorth(x);
        }

        // Fill in Y data
        for (int j = 0; j < N_cell; j++) {
            double y = grid.point(hps::Y, j);
            g0[j] = poisson.gEast(y);
            g1[j] = poisson.gWest(y);
            h0[j] = poisson.hEast(y);
            h1[j] = poisson.hWest(y);
        }

        // Stack boundary condition data
        g.intract(0*N_cell, g0);
        g.intract(1*N_cell, g1);
        g.intract(2*N_cell, g2);
        g.intract(3*N_cell, g3);
        h_exact.intract(0*N_cell, h0);
        h_exact.intract(1*N_cell, h1);
        h_exact.intract(2*N_cell, h2);
        h_exact.intract(3*N_cell, h3);

        // Fill in matrix data
        for (int i = 0; i < N_cell; i++) {
            for (int j = 0; j < N_cell; j++) {
                double x = grid.point(hps::X, i);
                double y = grid.point(hps::Y, j);
                u_exact.at(i,j) = poisson.u(x,y);
                f.at(i,j) = poisson.f(x,y);
            }
        }

        // Compute numerical solution to Poisson equation
        hps::Matrix<double> u_solve = mapSolution(grid, g, f);
        hps::Vector<double> h_solve = mapDirichletToNeumann(grid, g, f);

        // Compute errors and store
        double u_error = (u_solve - u_exact).gridNorm(grid.spacing[hps::X], grid.spacing[hps::Y], "infinity");
        double h_error = (h_solve - h_exact).gridNorm(grid.spacing[hps::X], grid.spacing[hps::Y], "infinity");
        u_errors[i] = u_error;
        h_errors[i] = h_error;

    }


}
