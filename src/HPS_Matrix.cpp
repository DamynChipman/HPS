#include "HPS_Matrix.hpp"

namespace hps {

/**
 * Multiplies a matrix `A` and a vector `x`. Wrapper for the BLAS dgemv.
 * @param   A  Matrix to multiply
 * @param   x  Vector to multiply
 * @return     New vector instance with result
 */
// hps::Vector<double> operator*(Matrix<double>& A, hps::Vector<double>& x) {
//
//     // Create output vector
//     hps::Vector<double> b(A.rows());
//
//     // Setup call
//     // void dgemv_(int* TRANS, int* M, int* N, double* ALPHA, double* A, int* LDA, double* X, int* INCX, double* BETA, double* Y, int* INCY);
//     char TRANS_ = 'C';
//     int M_ = A.rows();
//     int N_ = A.cols();
//     double ALPHA_ = 1.0;
//     double* A_ = A.data();
//     int LDA_ = A.cols();
//     double* X_ = x.data();
//     int INCX_ = 1;
//     double BETA_ = 0.0;
//     double* Y_ = b.data();
//     int INCY_ = 1;
//
//     dgemv_(&TRANS_, &M_, &N_, &ALPHA_, A_, &LDA_, X_, &INCX_, &BETA_, Y_, &INCY_);
//
//     return b;
// }


}
