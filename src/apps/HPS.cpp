#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include "HPS_Base.hpp"
#include "HPS_PatchSolver.hpp"
#include "HPS_Methods.hpp"
#include "HPS_IO.hpp"

using namespace hps;
using namespace std;

int main(int argc, char** argv) {

    // Parse inputs
    string usage = "HPS : <# of cells per side in leaf node> <poisson problem ID> <prefix for output files>";
    if (argc != 4) {
        cout << usage << endl;
        return EXIT_FAILURE;
    }

    int N = atoi(argv[1]);
    int prob_ID = atoi(argv[2]);
    string file_prefix = argv[3];

    // Create parent grid
    double lambda = 0.0;
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
    int N_levels = 3;
    int N_boundary = 2*N;
    N_pts[X] = 2*N; N_pts[Y] = 2*N;
    CellGrid<double, 2> domain_grid(N_pts, lower, upper);

    // Create Poisson problem
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

    // HPS Build stage
    cout << "================================== STARTING BUILD STAGE ==================================" << endl;
    int N_runs = 1;
    double build_time = 0;
    for (int n = 0; n < N_runs; n++) {
        cout << "----- run: " << n << " of " << N_runs << " -----" << endl;
        auto build_time_start = chrono::high_resolution_clock::now();
        buildStage(patches, N_levels, poisson, lambda);
        auto build_time_n = chrono::high_resolution_clock::now() - build_time_start;
        double build_time_seconds = (double) chrono::duration_cast<chrono::microseconds>(build_time_n).count() / 1e6;
        build_time = build_time + build_time_seconds;
    }
    build_time = build_time / N_runs;

    // HPS Solve stage
    cout << "================================== STARTING SOLVE STAGE ==================================" << endl;
    double solve_time = 0;
    for (int n = 0; n < N_runs; n++) {
        cout << "----- run: " << n << " of " << N_runs << " -----" << endl;
        auto solve_time_start = chrono::high_resolution_clock::now();
        solveStage(patches, N_levels, poisson, g_Dirichlet, lambda);
        auto solve_time_n = chrono::high_resolution_clock::now() - solve_time_start;
        double solve_time_seconds = (double) chrono::duration_cast<chrono::microseconds>(solve_time_n).count() / 1e6;
        solve_time = solve_time + solve_time_seconds;
    }
    solve_time = solve_time / N_runs;


    cout << "========================== STARTING REFINED PATCH SOLUTION STAGE =========================" << endl;
    Matrix<double> f_refined(2*N, 2*N);
    for (int i = 0; i < 2*N; i++) {
        for (int j = 0; j < 2*N; j++) {
            double x = refined_grid.point(X, i);
            double y = refined_grid.point(Y, j);
            f_refined.at(i,j) = poisson.f(x, y);
        }
    }
    Matrix<double> u_patch;
    double fishpack_time = 0;
    for (int n = 0; n < N_runs; n++) {
        cout << "----- run: " << n << " of " << N_runs << " -----" << endl;
        auto fishpack_time_start = chrono::high_resolution_clock::now();
        u_patch = mapSolution(refined_grid, g_Dirichlet, f_refined, lambda);
        auto fishpack_time_n = chrono::high_resolution_clock::now() - fishpack_time_start;
        double fishpack_time_seconds = (double) chrono::duration_cast<chrono::microseconds>(fishpack_time_n).count() / 1e6;
        fishpack_time = fishpack_time + fishpack_time_seconds;
    }
    fishpack_time = fishpack_time / N_runs;


    cout << "==================================== TIMING RESULTS ======================================" << endl;
    printf("Build Time: %12.6e seconds\n", build_time);
    printf("Solve Time: %12.6e seconds\n", solve_time);
    printf("Patch Time: %12.6e seconds\n", fishpack_time);

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
    // printf("Solution Inf-Norm Error: %8.4e\n", u_exact_error);
    printf("Patch Inf-Norm Error:    %8.4e\n", u_patch_error);
    // printf("FISHPACK Inf-Norm Error: %8.4e\n", u_fp_error);

    // Output results
    cout << "OUTPUTING DATA TO FILES WITH PREFIX: " << file_prefix << endl;
    string data_filename = file_prefix + "_data.txt";
    FILE* data_file = fopen(data_filename.c_str(), "w");
    fprintf(data_file, "Build Time: %12.6e seconds\n", build_time);
    fprintf(data_file, "Solve Time: %12.6e seconds\n", solve_time);
    fprintf(data_file, "Patch Time: %12.6e seconds\n", fishpack_time);
    fprintf(data_file, "# Cells [Leaf]:          %i\n", N);
    fprintf(data_file, "# Degrees of Freedom:    %i\n", (int) pow(2*N, 2));
    fprintf(data_file, "DtN Inf-Norm Error:      %8.4e\n", T_error);
    fprintf(data_file, "Solution Inf-Norm Error: %8.4e\n", u_exact_error);
    fprintf(data_file, "Patch Inf-Norm Error:    %8.4e\n", u_patch_error);
    fprintf(data_file, "FISHPACK Inf-Norm Error: %8.4e\n", u_fp_error);

    if (N <= 64) {
        string u_aprox_filename = file_prefix + "_u_aprox.vtk";
        string u_aprox_gridname = "refined_grid";
        string u_aprox_dataname = "u_aprox";
        writeVTK(refined_grid, u_aprox, u_aprox_filename, u_aprox_gridname, u_aprox_dataname);

        string u_exact_filename = file_prefix + "_u_exact.vtk";
        string u_exact_gridname = "refined_grid";
        string u_exact_dataname = "u_exact";
        writeVTK(refined_grid, u_exact, u_exact_filename, u_exact_gridname, u_exact_dataname);

        string u_patch_filename = file_prefix + "_u_patch.vtk";
        string u_patch_gridname = "refined_grid";
        string u_patch_dataname = "u_patch";
        writeVTK(refined_grid, u_patch, u_patch_filename, u_patch_gridname, u_patch_dataname);

        string u_exact_diff_filename = file_prefix + "_u_exact_diff.vtk";
        string u_exact_diff_gridname = "refined_grid";
        string u_exact_diff_dataname = "u_exact_diff";
        writeVTK(refined_grid, u_exact_diff, u_exact_diff_filename, u_exact_diff_gridname, u_exact_diff_dataname);

        string u_patch_diff_filename = file_prefix + "_u_patch_diff.vtk";
        string u_patch_diff_gridname = "refined_grid";
        string u_patch_diff_dataname = "u_patch_diff";
        writeVTK(refined_grid, u_patch_diff, u_patch_diff_filename, u_patch_diff_gridname, u_patch_diff_dataname);
    }

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

    return 0;
}
