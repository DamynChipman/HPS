#include "HPS_PatchSolver.hpp"

namespace hps {

hps::Matrix<double> mapSolution(hps::CellGrid<double, 2>& grid, hps::Vector<double>& g, hps::Matrix<double>& f, double lambda) {

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
    double ELMBDA = lambda;
    double* F = f_T.data();
    int IDIMF = M;
    double PERTRB;
    int IERROR;
    int WSIZE = 13*M + 4*N + M*((int)log2(N));
    double* W = (double*) malloc(WSIZE*sizeof(double));
    // std::cout << "WSIZE = " << WSIZE << std::endl;

    // std::cout << "A = " << A << std::endl;
    // std::cout << "B = " << B << std::endl;
    // std::cout << "M = " << M << std::endl;
    // std::cout << "C = " << C << std::endl;
    // std::cout << "D = " << D << std::endl;
    // std::cout << "N = " << N << std::endl;
    // std::cout << "MBDCND = " << MBDCND << std::endl;
    // std::cout << "NBDCND = " << NBDCND << std::endl;
    // std::cout << "ELMBDA = " << ELMBDA << std::endl;
    // for (int i = 0; i < N_side; i++) {
    //     printf("BDA[%i] = %8.4e    BDB[%i] = %8.4e    BDC[%i] = %8.4e    BDD[%i] = %8.4e\n", i, BDA[i], i, BDB[i], i, BDC[i], i, BDD[i]);
    // }
    // for (int i = 0; i < N_side; i++) {
    //     for (int j = 0; j < N_side; j++) {
    //         printf("F[%i][%i] = %8.4e\n", i, j, F[i*N_side + j]);
    //     }
    // }

    // Make Fortran call to FISHPACK
    // std::cout << "[mapSolution] beginning hstcrt_ call..." << std::endl;
    hstcrt_(&A, &B, &M, &MBDCND, BDA, BDB,
            &C, &D, &N, &NBDCND, BDC, BDD,
            &ELMBDA, F, &IDIMF, &PERTRB, &IERROR, W);
    // std::cout << "[mapSolution] end of hstcrt_ call..." << std::endl;
//     std::cout << "[mapSolution] IERROR = " << IERROR << std::endl;

    // std::cout << "A = " << A << std::endl;
    // std::cout << "B = " << B << std::endl;
    // std::cout << "M = " << M << std::endl;
    // std::cout << "C = " << C << std::endl;
    // std::cout << "D = " << D << std::endl;
    // std::cout << "N = " << N << std::endl;
    // std::cout << "MBDCND = " << MBDCND << std::endl;
    // std::cout << "NBDCND = " << NBDCND << std::endl;
    // std::cout << "ELMBDA = " << ELMBDA << std::endl;
    // for (int i = 0; i < N_side; i++) {
    //     printf("BDA[%i] = %8.4e    BDB[%i] = %8.4e    BDC[%i] = %8.4e    BDD[%i] = %8.4e\n", i, BDA[i], i, BDB[i], i, BDC[i], i, BDD[i]);
    // }
    // for (int i = 0; i < N_side; i++) {
    //     for (int j = 0; j < N_side; j++) {
    //         printf("F[%i][%i] = %8.4e\n", i, j, F[i*N_side + j]);
    //     }
    // }
    // std::cout << "here 1" << std::endl;
    if (IERROR != 0) {
        std::cerr << "[PatchSolver::mapSolution] WARNING: call to hstcrt_ returned non-zero error value: IERROR = " << IERROR << std::endl;
    }
    // std::cout << "IERROR = " << IERROR << std::endl;
    // std::cout << "here 2" << std::endl;
    // std::cout << "W[0] = " << W[0] << std::endl;

    // Move FISHPACK solution to matrix for output
    // std::cout << "here 3" << std::endl;
    hps::Matrix<double> u(grid.N_pts[hps::X], grid.N_pts[hps::Y]);
    // std::cout << "here 4" << std::endl;
    for (int i = 0; i < grid.N_pts[hps::X]; i++) {
        for (int j = 0; j < grid.N_pts[hps::Y]; j++) {
            // std::cout << "i = " << i << " j = " << j << std::endl;
            u.at(i,j) = F[j*grid.N_pts[hps::Y] + i];
        }
    }
    // std::cout << "here 5" << std::endl;

    return u;
}


hps::Vector<double> mapDirichletToNeumann(hps::CellGrid<double, 2> grid, hps::Vector<double> g, hps::Matrix<double> f, double lambda) {

    // Unpack Grid Data
    int N_cells = grid.N_pts[hps::X];

    // Unpack Dirichlet data
    hps::Vector<double> g0 = g.extract(0*N_cells, N_cells);
    hps::Vector<double> g1 = g.extract(1*N_cells, N_cells);
    hps::Vector<double> g2 = g.extract(2*N_cells, N_cells);
    hps::Vector<double> g3 = g.extract(3*N_cells, N_cells);

    // Compute Solution on Interior Nodes
    // std::cout << "mapping solution..." << std::endl;
    hps::Matrix<double> u = mapSolution(grid, g, f, lambda);
    // std::cout << "done with mapping solution..." << std::endl;
    // std::cout << "here now" << std::endl;

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


hps::Matrix<double> buildDirichletToNeumann(hps::CellGrid<double, 2> grid, double lambda) {

    // Variables
    int N = 4*grid.N_pts[hps::X];          // Number of rows/columns in T
    hps::Matrix<double> T(N, N);           // Output DtN map matrix T
    hps::Vector<double> e_hat_j(N, 0.0);   // Unit vector
    double z = 0.0;                        // Dummy variable for zero
    hps::Matrix<double> f_zero(N/4, N/4, z);   // Zero matrix for RHS function of PDE
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
        // std::cout << "here 1: " << j << std::endl;
        col_j = mapDirichletToNeumann(grid, e_hat_j, f_zero, lambda);

        // Intract j-th column into T
        // std::cout << "here 2: " << j << std::endl;
        T.intractColumn(j, col_j);
        // std::cout << "here " << j << std::endl;
    }

    return T;
}

} // Namespace hps
