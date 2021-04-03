#include <iostream>
#include "HPS_Base.hpp"
#include "HPS_Patch.hpp"
#include "HPS_PatchSolver.hpp"
#include "HPS_Merge.hpp"

using namespace std;
using namespace hps;

int main (int argc, char** argv) {

    double lambda = 0.0;

    int prob_ID;
    if (argc != 1) {
        prob_ID = atoi(argv[1]);
    }
    else {
        prob_ID = CONSTANT;
    }

    switch(prob_ID) {
		case CONSTANT:
			cout << "Solving uxx + uyy = 0\nu_exact = 1" << endl;
            break;
		case LINEAR:
			cout << "Solving uxx + uyy = 0\nu_exact = x + y" << endl;
            break;
		case QUAD:
			cout << "Solving uxx + uyy = 4\nu_exact = x^2 + y^2 + 2xy" << endl;
            break;
		case POLY:
			cout << "Solving uxx + uyy = 2x^2(6y^2 + x^2)\nu_exact = y^2 x^4" << endl;
            break;
		case TRIG:
			cout << "Solving uxx + uyy = -2(2 PI)^2 u\nu_exact = sin(2 PI x) sin(2 PI y)" << endl;
            break;
	}

    // Setup vectors with iteration info
    Vector<int> N_cells(5);
    int N_runs = N_cells.size();
    for (int i = 0; i < N_runs; i++) {
        N_cells[i] = pow(2, i+3);
    }
    Vector<double> T_errors(N_runs);
    Vector<double> h_errors_parent(N_runs);
    Vector<double> h_errors_merged(N_runs);

    cout << endl;
    cout << "N_PARENT: Number of cells in the parent grid (per side)" << endl;
    cout << "T_ERROR: Maximum difference between DtN matricies formed from parent and merged" << endl;
    cout << "H_ERROR_PARENT: Error between exact Neumann data and (T_parent * g + fhat)" << endl;
    cout << "H_ERROR_MERGED: Error between exact Neumann data and (T_merged * g + fhat)" << endl;
    cout << "========================= Error Analysis =========================" << endl;

    for (int i = 0; i < N_runs; i++) {

        // Setup problem inputs
        int N_cell = N_cells[i];
        int N_cell_leaf = N_cell / 2;
        int size = N_cell * N_cell;
        double Ax = -1.0;
        double Bx = 1.0;
        double Ay = -1.0;
        double By = 1.0;
        int Nx = N_cell;
        int Ny = N_cell;

        // Create finite difference grid
        int N_pts[2] = {Nx, Ny};
        double lower_bounds[2] = {Ax, Ay};
        double upper_bounds[2] = {Bx, By};
        CellGrid<double, 2> grid_parent(N_pts, lower_bounds, upper_bounds);

        // Create grids for patches
        //    Upper Left
        N_pts[X] = N_cell_leaf; N_pts[Y] = N_cell_leaf;
        lower_bounds[X] = Ax;   lower_bounds[Y] = (Ay + By)/2;
        upper_bounds[X] = (Ax + Bx)/2; upper_bounds[Y] = By;
        CellGrid<double, 2> grid_UL(N_pts, lower_bounds, upper_bounds);

        //    Upper Right
        lower_bounds[X] = (Ax + Bx)/2; lower_bounds[Y] = (Ay + By)/2;
        upper_bounds[X] = Bx;   upper_bounds[Y] = By;
        CellGrid<double, 2> grid_UR(N_pts, lower_bounds, upper_bounds);

        //    Lower Left
        lower_bounds[X] = Ax;   lower_bounds[Y] = Ay;
        upper_bounds[X] = (Ax + Bx)/2; upper_bounds[Y] = (Ay + By)/2;
        CellGrid<double, 2> grid_LL(N_pts, lower_bounds, upper_bounds);

        //    Lower Right
        lower_bounds[X] = (Ax + Bx)/2; lower_bounds[Y] = Ay;
        upper_bounds[X] = Bx;   upper_bounds[Y] = (Ay + By)/2;
        CellGrid<double, 2> grid_LR(N_pts, lower_bounds, upper_bounds);

        //    Merged upper
        N_pts[X] = 2*N_cell_leaf; N_pts[Y] = N_cell_leaf;
        lower_bounds[X] = Ax; lower_bounds[Y] = (Ay + By)/2;
        upper_bounds[X] = Bx; upper_bounds[Y] = By;
        CellGrid<double, 2> grid_upper(N_pts, lower_bounds, upper_bounds);

        //    Merged lower
        N_pts[X] = 2*N_cell_leaf; N_pts[Y] = N_cell_leaf;
        lower_bounds[X] = Ax; lower_bounds[Y] = Ay;
        upper_bounds[X] = Bx; upper_bounds[Y] = (Ay + By)/2;
        CellGrid<double, 2> grid_lower(N_pts, lower_bounds, upper_bounds);

        // Create Poisson problem
        PoissonProblem poisson(prob_ID, Ax, Bx, Ay, By);

        // Build grid_parent function vectors
        Vector<double> g0(N_cell);              // East Dirichlet data
        Vector<double> g1(N_cell);              // West Dirichlet data
        Vector<double> g2(N_cell);              // South Dirichlet data
        Vector<double> g3(N_cell);              // North Dirichlet data
        Vector<double> g(4*N_cell);             // Stacked Dirichlet data
        // Vector<double> g_lower(6*N_cell_leaf);
        // Vector<double> g_upper(6*N_cell_leaf);
        Vector<double> h0(N_cell);              // East Neumann data
        Vector<double> h1(N_cell);              // West Neumann data
        Vector<double> h2(N_cell);              // South Neumann data
        Vector<double> h3(N_cell);              // North Neumann data
        Vector<double> h_exact(4*N_cell);       // Stacked Neumann data
        Matrix<double> f(N_cell, N_cell);       // RHS matrix
        Matrix<double> f_UL(N_cell_leaf, N_cell_leaf);
        Matrix<double> f_UR(N_cell_leaf, N_cell_leaf);
        Matrix<double> f_LL(N_cell_leaf, N_cell_leaf);
        Matrix<double> f_LR(N_cell_leaf, N_cell_leaf);
        Matrix<double> u_exact(N_cell, N_cell); // Exact solution matrix

        // Fill in X data
        for (int i = 0; i < grid_parent.N_pts[X]; i++) {
            double x = grid_parent.point(X, i);
            g2[i] = poisson.u(x, Ay);
            g3[i] = poisson.u(x, By);
            h2[i] = poisson.dudy(x, Ay);
            h3[i] = poisson.dudy(x, By);
        }

        // Fill in Y data
        for (int j = 0; j < grid_parent.N_pts[Y]; j++) {
            double y = grid_parent.point(Y, j);
            g0[j] = poisson.u(Ax, y);
            g1[j] = poisson.u(Bx, y);
            h0[j] = poisson.dudx(Ax, y);
            h1[j] = poisson.dudx(Bx, y);
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
        for (int i = 0; i < grid_parent.N_pts[X]; i++) {
            for (int j = 0; j < grid_parent.N_pts[Y]; j++) {
                double x = grid_parent.point(X, i);
                double y = grid_parent.point(Y, j);
                u_exact.at(i,j) = poisson.u(x,y);
                f.at(i,j) = poisson.f(x,y);
            }
        }

        for (int i = 0; i < N_cell_leaf; i++) {
            for (int j = 0; j < N_cell_leaf; j++) {
                f_UL.at(i,j) = poisson.f(grid_UL.point(X, i), grid_UL.point(Y, j));
                f_UR.at(i,j) = poisson.f(grid_UR.point(X, i), grid_UR.point(Y, j));
                f_LL.at(i,j) = poisson.f(grid_LL.point(X, i), grid_LL.point(Y, j));
                f_LR.at(i,j) = poisson.f(grid_LR.point(X, i), grid_LR.point(Y, j));
            }
        }

        // Create level 0 T (parent)
        Matrix<double> T_parent = buildDirichletToNeumann(grid_parent, lambda);

        // Create patch data
        //    Matrix T
        Matrix<double> T_UL = buildDirichletToNeumann(grid_UL, lambda);
        Matrix<double> T_UR = buildDirichletToNeumann(grid_UR, lambda);
        Matrix<double> T_LL = buildDirichletToNeumann(grid_LL, lambda);
        Matrix<double> T_LR = buildDirichletToNeumann(grid_LR, lambda);

        //    Vector fhat
        Vector<double> g_zero_leaf(4*N_cell_leaf, 0.0);
        Vector<double> fhat_UL = mapDirichletToNeumann(grid_UL, g_zero_leaf, f_UL, lambda);
        Vector<double> fhat_UR = mapDirichletToNeumann(grid_UR, g_zero_leaf, f_UR, lambda);
        Vector<double> fhat_LL = mapDirichletToNeumann(grid_LL, g_zero_leaf, f_LL, lambda);
        Vector<double> fhat_LR = mapDirichletToNeumann(grid_LR, g_zero_leaf, f_LR, lambda);

        // if (N_cells[i] == 8) {
        //     cout << "f = \n" << f << endl;
        //     cout << "f_UL = \n" << f_UL << endl;
        //     cout << "f_UR = \n" << f_UR << endl;
        //     cout << "f_LL = \n" << f_LL << endl;
        //     cout << "f_LR = \n" << f_LR << endl;
        // }
        //
        // if (N_cells[i] == 8) {
        //     cout << "fhat_UL = \n" << fhat_UL << endl;
        //     cout << "fhat_UR = \n" << fhat_UR << endl;
        //     cout << "fhat_LL = \n" << fhat_LL << endl;
        //     cout << "fhat_LR = \n" << fhat_LR << endl;
        // }

        // Begin 4-to-1 merge from leaves
        //    Upper Horizontal Merge
        int N_levels = 3;
        int level = 2;
        merge_values upper_h_merge = mergeHorizontal(T_UL, T_UR, fhat_UL, fhat_UR, N_levels, level);
        Matrix<double> T_upper = upper_h_merge.T;
        Matrix<double> S_upper = upper_h_merge.S;
        Vector<double> fhat_upper = upper_h_merge.fhat;

        //    Lower Horizontal Merge
        merge_values lower_h_merge = mergeHorizontal(T_LL, T_LR, fhat_LL, fhat_LR, N_levels, level);
        Matrix<double> T_lower = lower_h_merge.T;
        Matrix<double> S_lower = lower_h_merge.S;
        Vector<double> fhat_lower = lower_h_merge.fhat;

        //    Vertical Merge
        level = 1;
        merge_values v_merge = mergeVertical(T_lower, T_upper, fhat_lower, fhat_upper, N_levels, level);
        Matrix<double> T_merged = v_merge.T;
        Matrix<double> S_merged = v_merge.S;
        Vector<double> fhat_merged = v_merge.fhat;

        if (N_cell == 8) {
            cout << "Size of S_merged: " << S_merged.rows() << " x " << S_merged.cols() << endl;
            cout << "Size of S_upper: " << S_upper.rows() << " x " << S_upper.cols() << endl;
            cout << "Size of S_lower: " << S_lower.rows() << " x " << S_lower.cols() << endl;
        }

        // // Create level 2 T (leaf)
        // int N_pts_leaf[2] = {N_cell_leaf, N_cell_leaf};
        // double lb_leaf[2] = {Ax, Ay};
        // double ub_leaf[2] = {Bx/2, By/2};
        // CellGrid<double, 2> grid_leaf(N_pts_leaf, lb_leaf, ub_leaf);
        // Matrix<double> T_level2 = buildDirichletToNeumann(grid_leaf, lambda);
        //
        // //    Horizontal merge
        // int N_levels = 3;
        // int level = 2;
        // merge_values h_merge = mergeHorizontal(T_level2, N_levels, level);
        // // auto [S_level1, T_level1] = mergeHorizontal(T_level2, N_levels, level);
        // // pair<Matrix<double>, Matrix<double>> ST_pair_level1 = mergeHorizontal(T_level2, N_levels, level);
        // Matrix<double> S_level1 = h_merge.S;
        // Matrix<double> T_level1 = h_merge.T;
        //
        // //    Vertical merge
        // level = 1;
        // merge_values v_merge = mergeVertical(T_level1, N_levels, level);
        // // auto [S_level0, T_level0_merged] = mergeVertical(T_level1, N_levels, level);
        // // pair<Matrix<double>, Matrix<double>> ST_pair_level0 = mergeVertical(T_level1, N_levels, level);
        // Matrix<double> S_level0 = v_merge.S;
        // Matrix<double> T_level0_merged = v_merge.T;

        // Verifications
        //    Entry differences
        double T_error = 0.0;
        for (int i = 0; i < 4*N_cell; i++) {
            for (int j = 0; j < 4*N_cell; j++) {
                T_error = max(T_error, fabs(T_merged.at(i,j) - T_parent.at(i,j)));
            }
        }
        T_errors[i] = T_error;

        //     Mapping
        Vector<double> g_zero(4*N_cell, 0.0);
        Vector<double> fhat = mapDirichletToNeumann(grid_parent, g_zero, f, lambda); // Does T * g -> h


        // Tg + fhat = h



        Vector<double> h_mapped_parent = T_parent * g + fhat;
        Vector<double> h_mapped_merged = T_merged * g + fhat_merged;
        double h_error_parent = 0.0;
        double h_error_merged = 0.0;
        double h_error_diff = 0.0;
        double fhat_error_diff = 0.0;
        for (int i = 0; i < 4*N_cell; i++) {
            h_error_parent = max(h_error_parent, fabs(h_mapped_parent[i] - h_exact[i]));
            h_error_merged = max(h_error_merged, fabs(h_mapped_merged[i] - h_exact[i]));
            h_error_diff = max(h_error_diff, fabs(h_mapped_parent[i] - h_mapped_merged[i]));
            fhat_error_diff = max(fhat_error_diff, fabs(fhat[i] - fhat_merged[i]));
        }
        h_errors_parent[i] = h_error_parent;
        h_errors_merged[i] = h_error_merged;
        // cout << "FHAT_ERROR_DIFF = " << fhat_error_diff << endl;

        // if (N_cells[i] == 8) {
        //     cout << "fhat = \n" << fhat << endl;
        //     cout << "fhat_merged = \n" << fhat_merged << endl;
        //     cout << "h_mapped_parent = \n" << h_mapped_parent << endl;
        //     cout << "h_mapped_merged = \n" << h_mapped_merged << endl;
        //     cout << "h_diff = \n" << (h_mapped_merged - h_mapped_parent).abs() << endl;
        // }

        // Output results
        printf("N_PARENT = %6i    T_ERROR = %8.4e    H_ERROR_PARENT = %8.4e    H_ERROR_MERGED = %8.4e    H_ERROR_DIFF = %8.4e\n", N_cells[i], T_errors[i], h_errors_parent[i], h_errors_merged[i], h_error_diff);

    }

    return 0;
}
