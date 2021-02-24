#include "HPS_Base.hpp"
#include "HPS_PatchSolver.hpp"

#include <iostream>
#include <cmath>

using namespace hps;
using namespace std;

int main(int argc, char* argv[]) {

    double lambda = 0.0; // Zero for Poisson's Equation (and not Helmholtz)

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
    Vector<int> N_cells(7);
    // Vector<int> N_cells({2048});
    // Vector<int> N_cells({100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000});
    // Vector<int> N_cells({1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000});
    // Vector<int> N_cells({1500, 1510, 1520, 1530, 1540, 1550, 1560, 1570, 1580, 1590, 1600});
    // Vector<int> N_cells({1530, 1531, 1532, 1533, 1534, 1535, 1536, 1537, 1538, 1539, 1540});
    int N_runs = N_cells.size();
    Vector<double> u_errors(N_runs);
    Vector<double> h_errors(N_runs);
    Vector<double> T_errors(N_runs);
    Vector<double> delta_xs(N_runs);
    Vector<double> delta_ys(N_runs);
    for (int i = 0; i < N_runs; i++) {
        N_cells[i] = pow(2, i+2);
    }

    // Begin convergence analysis loop
    for (int i = 0; i < N_runs; i++) {

        // Setup problem inputs
        int N_cell = N_cells[i];
        // cout << "Running tests for N = " << N_cell << endl;
        int size = N_cell * N_cell;
        double Ax = 0.0;
        double Bx = 1.0;
        double Ay = 0.0;
        double By = 1.0;
        int Nx = N_cell;
        int Ny = N_cell;

        // Create finite difference grid
        int N_pts[2] = {Nx, Ny};
        double lower_bounds[2] = {Ax, Ay};
        double upper_bounds[2] = {Bx, By};
        CellGrid<double, 2> grid(N_pts, lower_bounds, upper_bounds);

        // Create Poisson problem
        PoissonProblem poisson(prob_ID, Ax, Bx, Ay, By);

        // Build grid function vectors
        Vector<double> g0(N_cell);              // East Dirichlet data
        Vector<double> g1(N_cell);              // West Dirichlet data
        Vector<double> g2(N_cell);              // South Dirichlet data
        Vector<double> g3(N_cell);              // North Dirichlet data
        Vector<double> g(4*N_cell);             // Stacked Dirichlet data
        Vector<double> h0(N_cell);              // East Neumann data
        Vector<double> h1(N_cell);              // West Neumann data
        Vector<double> h2(N_cell);              // South Neumann data
        Vector<double> h3(N_cell);              // North Neumann data
        Vector<double> h_exact(4*N_cell);       // Stacked Neumann data
        Matrix<double> f(N_cell, N_cell);       // RHS matrix
        Matrix<double> u_exact(N_cell, N_cell); // Exact solution matrix

        // Fill in X data
        for (int i = 0; i < grid.N_pts[X]; i++) {
            double x = grid.point(X, i);
            g2[i] = poisson.u(x, Ay);
            g3[i] = poisson.u(x, By);
            h2[i] = poisson.dudy(x, Ay);
            h3[i] = poisson.dudy(x, By);
        }

        // Fill in Y data
        for (int j = 0; j < grid.N_pts[Y]; j++) {
            double y = grid.point(Y, j);
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
        for (int i = 0; i < grid.N_pts[X]; i++) {
            for (int j = 0; j < grid.N_pts[Y]; j++) {
                double x = grid.point(X, i);
                double y = grid.point(Y, j);
                u_exact.at(i,j) = poisson.u(x,y);
                f.at(i,j) = poisson.f(x,y);
            }
        }

        // Compute numerical solution to Poisson equation
        Matrix<double> u_solve = mapSolution(grid, g, f, lambda);
        Vector<double> h_solve = mapDirichletToNeumann(grid, g, f, lambda);

        // Build DtN matrix
        Matrix<double> T = buildDirichletToNeumann(grid, lambda);
        Vector<double> g_zero(4*N_cell, 0.0);
        Vector<double> f_hat = mapDirichletToNeumann(grid, g_zero, f, lambda);
        Vector<double> h_mapped = T * g + f_hat;

        // Compute errors and store
        double delta_x = grid.spacing[X];
        double delta_y = grid.spacing[Y];
        double u_error = 0.0;
        double h_error = 0.0;
        double T_error = 0.0;
        for (int i = 0; i < N_cell; i++) {
            for (int j = 0; j < N_cell; j++) {
                u_error = max(u_error, fabs(u_solve.at(i,j) - u_exact.at(i,j)));
            }
        }
        for (int i = 0; i < 4*N_cell; i++) {
            h_error = max(h_error, fabs(h_solve[i] - h_exact[i]));
        }
        for (int i = 0; i < 4*N_cell; i++) {
            T_error = max(T_error, fabs(h_mapped[i] - h_exact[i]));
        }
        // double u_error = (u_solve - u_exact).max();
        // double h_error = (h_solve - h_exact).max();
        // double u_error = (u_exact - u_solve).gridNorm(delta_x, delta_y, 2);
        // double h_error = (h_exact - h_solve).gridNorm(delta_x, delta_y, 2);
        u_errors[i] = u_error;
        h_errors[i] = h_error;
        T_errors[i] = T_error;
        delta_xs[i] = delta_x;
        delta_ys[i] = delta_y;

        // Print errors
        printf("N = %6i    U_ERROR = %16.8e    H_ERROR = %16.8e    T_ERROR = %16.8e\n", N_cells[i], u_errors[i], h_errors[i], T_errors[i]);
    }
    // cout << "Done! Printing results..." << endl;

    // Output results
    // for (int i = 0; i < N_runs; i++) {
    //     printf("N = %6i    U_ERROR = %16.8e    H_ERROR = %16.8e    T_ERROR = %16.8e\n", N_cells[i], u_errors[i], h_errors[i], T_errors[i]);
    // }

    for (int i = 0; i < N_runs - 1; i++) {
        double ux_order = (log(u_errors[i]/u_errors[i+1])) / (log(delta_xs[i]/delta_xs[i+1]));
        double uy_order = (log(u_errors[i]/u_errors[i+1])) / (log(delta_ys[i]/delta_ys[i+1]));
        double hx_order = (log(h_errors[i]/h_errors[i+1])) / (log(delta_xs[i]/delta_xs[i+1]));
        double hy_order = (log(h_errors[i]/h_errors[i+1])) / (log(delta_ys[i]/delta_ys[i+1]));
        printf("ux_order = %16.8f    uy_order = %16.8f    hx_order = %16.8f    hy_order = %16.8f\n", ux_order, uy_order, hx_order, hy_order);
    }

    // Why does the problem diverge from the solution past N_cells > 1536???

    return 0;

}
