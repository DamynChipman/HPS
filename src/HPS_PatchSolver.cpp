#include "HPS_PatchSolver.hpp"

namespace hps {

hps::Matrix<double> mapSolution(hps::CellGrid<double, 2>& grid, hps::Vector<double>& g, hps::Matrix<double>& f) {

    // Unpack Dirichlet data
//     std::cout << "[mapSolution] unpacking Dirichlet data..." << std::endl;
    int N_side = grid.N_pts[hps::X];
    hps::Vector<double> g_Ax = g.extract(0*N_side, N_side);
    hps::Vector<double> g_Bx = g.extract(1*N_side, N_side);
    hps::Vector<double> g_Ay = g.extract(2*N_side, N_side);
    hps::Vector<double> g_By = g.extract(3*N_side, N_side);

    // Transpose RHS for Fortran call
    hps::Matrix<double> f_T = f.transpose();

    // Setup Fortran call to FISHPACK
//     std::cout << "[mapSolution] setting up Fortran call..." << std::endl;
    double A = grid.lower_limit[hps::X];
    double B = grid.upper_limit[hps::X];
    int M = grid.N_pts[hps::X];
    int MBDCND = 1;
    double* BDA = g_Ax.data();
    double* BDB = g_Bx.data();
    double C = grid.lower_limit[hps::Y];
    double D = grid.upper_limit[hps::Y];
    int N = grid.N_pts[hps::Y];
    int NBDCND = 1;
    double* BDC = g_Ay.data();
    double* BDD = g_By.data();
    double ELMBDA = 0.0;
    double* F = f_T.data();
    int IDIMF = M;
    double PERTRB;
    int IERROR;

    // Make Fortran call to FISHPACK
//     std::cout << "[mapSolution] beginning hstcrt_ call..." << std::endl;
    hstcrt_(&A, &B, &M, &MBDCND, BDA, BDB,
            &C, &D, &N, &NBDCND, BDC, BDD,
            &ELMBDA, F, &IDIMF, &PERTRB, &IERROR);
//     std::cout << "[mapSolution] end of hstcrt_ call..." << std::endl;
//     std::cout << "[mapSolution] IERROR = " << IERROR << std::endl;

    if (IERROR != 0) {
        std::cerr << "[PatchSolver::mapSolution] WARNING: call to hstcrt_ returned non-zero error value: IERROR = " << IERROR << std::endl;
    }

    // Move FISHPACK solution to matrix for output
    hps::Matrix<double> u(f.rows(), f.cols());
    for (int i = 0; i < grid.N_pts[hps::X]; i++) {
        for (int j = 0; j < grid.N_pts[hps::Y]; j++) {
            u.at(i,j) = F[j*grid.N_pts[hps::Y] + i];
        }
    }

    return u;
}


hps::Vector<double> mapDirichletToNeumann(hps::CellGrid<double, 2> grid, hps::Vector<double> g, hps::Matrix<double> f) {

    // Unpack Grid Data
    int N_cells = grid.N_pts[hps::X];

    // Unpack Dirichlet data
    hps::Vector<double> g0 = g.extract(0*N_cells, N_cells);
    hps::Vector<double> g1 = g.extract(1*N_cells, N_cells);
    hps::Vector<double> g2 = g.extract(2*N_cells, N_cells);
    hps::Vector<double> g3 = g.extract(3*N_cells, N_cells);

    // Compute Solution on Interior Nodes
    hps::Matrix<double> u = mapSolution(grid, g, f);

    // Get Interior Edge Cell Data and Compute Neumann Data
    //    Interior Cell Data
    hps::Vector<double> u0 = u.extractRow(0);
    hps::Vector<double> u1 = u.extractRow(N_cells-1);
    hps::Vector<double> u2 = u.extractColumn(0);
    hps::Vector<double> u3 = u.extractColumn(N_cells-1);

    //    Neumann Data
    double dtn_x = 2.0 / grid.spacing[hps::X];
    double dtn_y = 2.0 / grid.spacing[hps::Y];
    hps::Vector<double> h0 = (u0 - g0) * (dtn_x);
    hps::Vector<double> h1 = (u1 - g1) * (-dtn_x);
    hps::Vector<double> h2 = (u2 - g2) * (dtn_y);
    hps::Vector<double> h3 = (u3 - g3) * (-dtn_y);

    // Column stack Neumann data
    hps::Vector<double> h(4*N_cells);
    h.intract(0*N_cells, h0);
    h.intract(1*N_cells, h1);
    h.intract(2*N_cells, h2);
    h.intract(3*N_cells, h3);

    // if (N_cells == 8) {
    //     std::cout << "u = " << u << std::endl;
    //     std::cout << "u0 = " << u0 << std::endl;
    //     std::cout << "u1 = " << u1 << std::endl;
    //     std::cout << "u2 = " << u2 << std::endl;
    //     std::cout << "u3 = " << u3 << std::endl;
    //     std::cout << "g0 = " << g0 << std::endl;
    //     std::cout << "g1 = " << g1 << std::endl;
    //     std::cout << "g2 = " << g2 << std::endl;
    //     std::cout << "g3 = " << g3 << std::endl;
    //     std::cout << "h0 = " << h0 << std::endl;
    //     std::cout << "h1 = " << h1 << std::endl;
    //     std::cout << "h2 = " << h2 << std::endl;
    //     std::cout << "h3 = " << h3 << std::endl;
    // }

    return h;
}


hps::Matrix<double> buildDirichletToNeumann(hps::CellGrid<double, 2> grid) {

    // Variables
    int N = 4*grid.N_pts[hps::X];          // Number of rows/columns in T
    hps::Matrix<double> T(N, N);           // Output DtN map matrix T
    hps::Vector<double> e_hat_j(N, 0.0);   // Unit vector
    double z = 0.0;                        // Dummy variable for zero
    hps::Matrix<double> f_zero(N, N, z);   // Zero matrix for RHS function of PDE
    hps::Vector<double> col_j(N);          // j-th column of T

    // Iterate through columns of identity matrix with e_hat_j
    for (int j = 0; j < N; j++) {

        // Correct unit vector
        if (j == 0) {
            e_hat_j[j] = 1.0;
        }
        else {
            e_hat_j[j-1] = 0.0;
            e_hat_j[j] = 1.0;
        }

        // Perform DtN map on unit vector to get i-th column of T
        col_j = mapDirichletToNeumann(grid, e_hat_j, f_zero);

        // Intract j-th column into T
        T.intractColumn(j, col_j);
    }

    return T;
}

} // Namespace hps
