#include <iostream>
#include <chrono>
#include <ctime>
#include "HPS_Base.hpp"
#include "HPS_PatchSolver.hpp"
// #include "HPS_Merge.hpp"
#include "HPS_Methods.hpp"

using namespace hps;
using namespace std;

int main(int argc, char** argv) {

    // Debugging
    // cout << "----- debugging area -----" << endl;
    // int r1 = 2;
    // int c1 = 8;
    // int r2 = 8;
    // int c2 = 8;
    //
    // hps::Matrix<double> S(r1, c1,
    //     {
    //         1, 2, 3, 4, 5, 6, 7, 8,
    //         10, 20, 30, 40, 50, 60, 70, 80
    //     }
    // );
    // hps::Matrix<double> P(r2, c2,
    //     {
    //         1, 0, 0, 0, 0, 0, 0, 0,
    //         0, 1, 0, 0, 0, 0, 0, 0,
    //         0, 0, 0, 0, 1, 0, 0, 0,
    //         0, 0, 0, 0, 0, 1, 0, 0,
    //         0, 0, 1, 0, 0, 0, 0, 0,
    //         0, 0, 0, 1, 0, 0, 0, 0,
    //         0, 0, 0, 0, 0, 0, 1, 0,
    //         0, 0, 0, 0, 0, 0, 0, 1,
    //     }
    // );
    //
    // cout << "S = \n" << S << endl;
    // cout << "P = \n" << P << endl;
    //
    // hps::Matrix<double> S_test = S * P;
    //
    // cout << "S_test = \n" << S_test << endl;
    //
    // cout << "----- END OF DEBUGGING -----" << endl;
    // return 0;

    // Create parent grid
    double lambda = 0;
    int N = atoi(argv[1]);
    double xL = -1;
    double xU = 1;
    double yL = -1;
    double yU = 1;
    int N_pts[2] = {N, N};
    double lower[2] = {xL, yL};
    double upper[2] = {xU, yU};
    CellGrid<double, 2> parent_grid(N_pts, lower, upper);

    // Create global grid
    int N_refined[2] = {2*N, 2*N};
    CellGrid<double, 2> refined_grid(N_refined, lower, upper);
    Patch parent_patch(refined_grid);
    parent_patch.buildT(lambda);

    // Create domain grid
    int N_levels = atoi(argv[2]);
    int N_boundary = 2*N;
    N_pts[X] = 2*N; N_pts[Y] = 2*N;
    CellGrid<double, 2> domain_grid(N_pts, lower, upper);

    // Create Poisson problem
    int prob_ID = atoi(argv[3]);
    PoissonProblem poisson(prob_ID, xL, xU, yL, yU);

    //    Dirichlet boundary data
    Vector<double> g0(N_boundary);
    Vector<double> g1(N_boundary);
    Vector<double> g2(N_boundary);
    Vector<double> g3(N_boundary);
    Vector<double> g_Dirichlet(4*N_boundary);

    //    Fill in X data
    for (int i = 0; i < domain_grid.N_pts[X]; i++) {
        double x = domain_grid.point(X, i);
        g2[i] = poisson.u(x, yL);
        g3[i] = poisson.u(x, yU);
    }

    //    Fill in Y data
    for (int j = 0; j < domain_grid.N_pts[Y]; j++) {
        double y = domain_grid.point(Y, j);
        g0[j] = poisson.u(xL, y);
        g1[j] = poisson.u(xU, y);
    }

    //    Stack boundary data
    g_Dirichlet.intract(0*N_boundary, g0);
    g_Dirichlet.intract(1*N_boundary, g1);
    g_Dirichlet.intract(2*N_boundary, g2);
    g_Dirichlet.intract(3*N_boundary, g3);

    //    Fill in matrix data
    Matrix<double> u_exact(N_boundary, N_boundary);
    for (int i = 0; i < domain_grid.N_pts[X]; i++) {
        for (int j = 0; j < domain_grid.N_pts[Y]; j++) {
            double x = domain_grid.point(X, i);
            double y = domain_grid.point(Y, j);
            u_exact.at(i,j) = poisson.u(x,y);
        }
    }

    // Setup HPS
    int N_patches = computeNumberOfPatches(N_levels);
    Vector<Patch> patches = setupHPS(parent_grid, N_levels);

    for (int tau = 0; tau < N_patches; tau++) {
        cout << "----------------------------" << endl;
        cout << "tau = " << tau << endl;
        cout << "ID = " << patches[tau].ID << endl;
        cout << "level = " << patches[tau].level << endl;
        cout << "coords = [" << patches[tau].coords[0] << ", " << patches[tau].coords[1] << "]" << endl;
        cout << "grid = " << endl;
        cout << "    Nx = " << patches[tau].grid.N_pts[X] << endl;
        cout << "    Ny = " << patches[tau].grid.N_pts[Y] << endl;
        cout << "    xL = " << patches[tau].grid.lower_limit[X] << endl;
        cout << "    xU = " << patches[tau].grid.upper_limit[X] << endl;
        cout << "    yL = " << patches[tau].grid.lower_limit[Y] << endl;
        cout << "    yU = " << patches[tau].grid.upper_limit[Y] << endl;
        cout << endl;
    }

    // for (int tau = N_patches-1; tau >= 0; tau--) {
    //     cout << tau << "    " << patches[tau].ID << "    " << patches[tau].level << "    " << patches[tau].is_leaf << endl;
    // }

    // HPS Build stage
    cout << "================================== STARTING BUILD STAGE ==================================" << endl;
    clock_t build_time;
    build_time = clock();
    buildStage(patches, N_levels, poisson, lambda);
    build_time = clock() - build_time;

    // HPS Solve stage
    cout << "================================== STARTING SOLVE STAGE ==================================" << endl;
    clock_t solve_time;
    solve_time = clock();
    solveStage(patches, N_levels, poisson, g_Dirichlet, lambda);
    solve_time = clock() - solve_time;
    // for (int i = 0; i < domain_grid.N_pts[X]; i++) {
    //     cout << poisson.u(domain_grid.point(X,i), 0) << ",  ";
    // }
    // cout << endl;

    // cout << "========================== STARTING REFINED PATCH SOLUTION STAGE =========================" << endl;
    Matrix<double> f_refined(2*N, 2*N);
    for (int i = 0; i < 2*N; i++) {
        for (int j = 0; j < 2*N; j++) {
            double x = refined_grid.point(X, i);
            double y = refined_grid.point(Y, j);
            f_refined.at(i,j) = poisson.f(x, y);
        }
    }
    clock_t fishpack_time = clock();
    Matrix<double> u_patch = mapSolution(refined_grid, g_Dirichlet, f_refined, lambda);
    fishpack_time = clock() - fishpack_time;

    cout << "==================================== TIMING RESULTS ======================================" << endl;
    printf("Build Time: %f seconds\n", ((float)build_time) / CLOCKS_PER_SEC);
    printf("Solve Time: %f seconds\n", ((float)solve_time) / CLOCKS_PER_SEC);
    printf("Patch Time: %f seconds\n", ((float)fishpack_time) / CLOCKS_PER_SEC);

    // Check error


    // DtN Error
    Matrix<double> T_merged = patches[0].T;
    Matrix<double> T_parent = parent_patch.T;
    Matrix<double> T_diff = T_merged - T_parent;
    double T_error = 0;
    for (int i = 0; i < T_merged.rows(); i++) {
        for (int j = 0; j < T_merged.cols(); j++) {
            T_error = max(T_error, fabs(T_diff.at(i,j)));
        }
    }
    // test

    int N_half_boundary = N_boundary / 2;
    Matrix<double> u_aprox(N_boundary, N_boundary);
    u_aprox.intract(0*N_half_boundary, 0*N_half_boundary, patches[3].u);
    u_aprox.intract(0*N_half_boundary, 1*N_half_boundary, patches[5].u);
    u_aprox.intract(1*N_half_boundary, 0*N_half_boundary, patches[4].u);
    u_aprox.intract(1*N_half_boundary, 1*N_half_boundary, patches[6].u);

    Matrix<double> u_exact_diff = u_exact - u_aprox;
    Matrix<double> u_patch_diff = u_patch - u_aprox;
    Matrix<double> u_fp_diff = u_exact - u_patch;

    double u_exact_error = 0.0;
    double u_patch_error = 0.0;
    double u_fp_error = 0.0;
    for (int i = 0; i < N_boundary; i++) {
        for (int j = 0; j < N_boundary; j++) {
            u_exact_error = max(u_exact_error, fabs(u_aprox.at(i,j) - u_exact.at(i,j)));
            u_patch_error = max(u_patch_error, fabs(u_aprox.at(i,j) - u_patch.at(i,j)));
            u_fp_error = max(u_fp_error, fabs(u_exact.at(i,j) - u_patch.at(i,j)));
        }
    }

    cout << "===================================== ERROR RESULTS ======================================" << endl;
    switch(prob_ID) {
		case CONSTANT:
			cout << "Solving uxx + uyy = 0\nu_exact = 1" << endl;
            break;
		case LINEAR:
			cout << "Solving uxx + uyy = 0\nu_exact = x + y" << endl;
            break;
        case LAPLACE:
            cout << "Solving uxx + uyy = 0\nu_exact = y*sin(2*M_PI*x) + x*cos(2*M_PI*y) + 4" << endl;
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
    printf("# Cells [Leaf]:          %i\n", N);
    printf("# Degrees of Freedom:    %i\n", (int) pow(2*N, 2));
    printf("DtN Inf-Norm Error:      %8.4e\n", T_error);
    printf("Solution Inf-Norm Error: %8.4e\n", u_exact_error);
    printf("Patch Inf-Norm Error:    %8.4e\n", u_patch_error);
    printf("FISHPACK Inf-Norm Error: %8.4e\n", u_fp_error);
    if (N < 10) {
        cout << "u_exact = \n" << u_exact << endl;
        cout << "u_patch = \n" << u_patch << endl;
        cout << "u_aprox = \n" << u_aprox << endl;
        cout << "u_exact_diff = \n" << u_exact_diff << endl;

        cout << "----- Level 0 -> Level 1 -----" << endl;
        cout << "patches[0].g = \n" << patches[0].g << endl;
        cout << "patches[1].g = \n" << patches[1].g << endl;
        cout << "patches[2].g = \n" << patches[2].g << endl;

        Vector<double> g_interior(N_boundary);
        for (int i = 0; i < domain_grid.N_pts[X]; i++) {
            double x = domain_grid.point(X, i);
            double y = 0.0;
            g_interior[i] = poisson.u(x, y);
        }
        Vector<double> g_mapped = patches[0].S * g_Dirichlet;
        cout << "g_interior = \n" << g_interior << endl;
        cout << "g_mapped = \n" << g_mapped << endl;

        cout << "----- Level 1 -> Level 2 -----" << endl;
        cout << "patches[1].g = \n" << patches[1].g << endl;
        cout << "patches[3].g = \n" << patches[3].g << endl;
        cout << "patches[4].g = \n" << patches[4].g << endl;

        g_interior = Vector<double>(N);
        for (int j = 0; j < domain_grid.N_pts[Y]/2; j++) {
            double x = 0;
            double y = domain_grid.point(Y, j);
            g_interior[j] = poisson.u(x, y);
        }
        g_mapped = patches[1].S * patches[1].g;
        cout << "g_interior = \n" << g_interior << endl;
        cout << "g_mapped = \n" << g_mapped << endl;

    }


    // cout << "S = \n" << patches[0].S << endl;

    // cout << "patches[3].u = " << endl << patches[3].u << endl;
    // cout << "patches[4].u = " << endl << patches[4].u << endl;
    // cout << "patches[5].u = " << endl << patches[5].u << endl;
    // cout << "patches[6].u = " << endl << patches[6].u << endl;
    //
    // cout << "===== ===== =====" << endl;
    // cout << "u_exact = " << endl << u_exact << endl;


    // Matrix<double> T_merged = patches[0].T;
    // Matrix<double> T_parent =


    return 0;
}
