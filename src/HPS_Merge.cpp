#include "HPS_Merge.hpp"

namespace hps {

/**
 * Performs the merge operation (linear algebra) for the solution operator S.
 * @NOTE: This HAS been unit tested against Python and passed
 * @param  T_alpha_33 Matrix for T_alpha_33
 * @param  T_beta_33  Matrix for T_beta_33
 * @param  T_alpha_31 Matrix for T_alpha_31
 * @param  T_beta_32  Matrix for T_beta_32
 * @return            Computed solution operator S
 */
Matrix<double> mergeOperationS(
    Matrix<double> T_alpha_33,
    Matrix<double> T_beta_33,
    Matrix<double> T_alpha_31,
    Matrix<double> T_beta_32) {

        Matrix<double> S_RHS(T_alpha_31.rows(), T_alpha_31.cols() + T_beta_32.cols());
        T_alpha_31.negate();
        S_RHS.intract(0, 0, T_alpha_31);
        S_RHS.intract(0, T_alpha_31.cols(), T_beta_32);
        Matrix<double> S = solve(T_alpha_33 - T_beta_33, S_RHS);
        return S;

}

/**
 * Performs the merge operation (linear algebra) for the DtN operator T.
 * @NOTE: This HAS been unit tested against Python and passed
 * @param  T_alpha_11 Matrix for T_alpha_11
 * @param  T_beta_22  Matrix for T_beta_22
 * @param  T_alpha_13 Matrix for T_alpha_13
 * @param  T_beta_23  Matrix for T_beta_23
 * @param  S          Matrix for solution operator
 * @return            Computed DtN operator T
 */
Matrix<double> mergeOperationT(
    Matrix<double> T_alpha_11,
    Matrix<double> T_beta_22,
    Matrix<double> T_alpha_13,
    Matrix<double> T_beta_23,
    Matrix<double> S) {

        double z = 0.0;
        Matrix<double> T(T_alpha_11.rows() + T_beta_22.rows(), T_alpha_11.cols() + T_beta_22.cols(), z);
        T.intract(0, 0, T_alpha_11);
        T.intract(T_alpha_11.rows(), T_alpha_11.cols(), T_beta_22);
        Matrix<double> T_RHS(T_alpha_13.rows() + T_beta_23.rows(), T_alpha_13.cols());
        T_RHS.intract(0, 0, T_alpha_13);
        T_RHS.intract(T_alpha_13.rows(), 0, T_beta_23);
        T += T_RHS*S;
        return T;

}

/**
 * Performs the merge operation (linear algebra) for the body load vector fhat
 * @param  T_alpha_13   Matrix for T_alpha_13
 * @param  T_beta_23    Matrix for T_beta_23
 * @param  T_alpha_33   Matrix for T_alpha_33
 * @param  T_beta_33    Matrix for T_beta_33
 * @param  fhat_alpha_1 Vector for fhat_alpha_1
 * @param  fhat_beta_2  Vector for fhat_beta_2
 * @param  fhat_alpha_3 Vector for fhat_alpha_3
 * @param  fhat_beta_3  Vector for fhat_beta_3
 * @return              Computed fhat vector
 */
Vector<double> mergeOperationF(
    Matrix<double> T_alpha_13,
    Matrix<double> T_beta_23,
    Matrix<double> T_alpha_33,
    Matrix<double> T_beta_33,
    Vector<double> fhat_alpha_1,
    Vector<double> fhat_beta_2,
    Vector<double> fhat_alpha_3,
    Vector<double> fhat_beta_3) {

        Vector<double> fhat(fhat_alpha_1.size() + fhat_beta_2.size());
        fhat.intract(0, fhat_alpha_1);
        fhat.intract(fhat_alpha_1.size(), fhat_beta_2);
        Matrix<double> T_RHS(T_alpha_13.rows() + T_beta_23.rows(), T_alpha_13.cols());
        T_RHS.intract(0, 0, T_alpha_13);
        T_RHS.intract(T_alpha_13.rows(), 0, T_beta_23);
        Vector<double> T_solved = solve(T_alpha_33 - T_beta_33, fhat_beta_3 - fhat_alpha_3);
        fhat += T_RHS * T_solved;
        return fhat;

}

/**
 * Performs the merge operation (linear algebra) for the w vector
 * @param  T_alpha_33   Matrix for T_alpha_33
 * @param  T_beta_33    Matrix for T_beta_33
 * @param  fhat_alpha_3 Vector for fhat_alpha_3
 * @param  fhat_beta_3  Vector for fhat_beta_3
 * @return              Computed w vector
 */
Vector<double> mergeOperationW(
    Matrix<double> T_alpha_33,
    Matrix<double> T_beta_33,
    Vector<double> fhat_alpha_3,
    Vector<double> fhat_beta_3) {

        Vector<double> w = solve(T_alpha_33 - T_beta_33, fhat_beta_3 - fhat_alpha_3);
        return w;

}

/**
 * Performs a horizontal merge.
 * @param  T_alpha    Alpha DtN matrix
 * @param  T_beta     Beta DtN matrix
 * @param  fhat_alpha Alpha fhat vector
 * @param  fhat_beta  Beta fhat vector
 * @param  N_levels   Number of levels in tree
 * @param  level      Current level in tree
 * @return            merge_values struct with parent patch data
 */
merge_values mergeHorizontal(
    Matrix<double>& T_alpha,
    Matrix<double>& T_beta,
    Vector<double>& fhat_alpha,
    Vector<double>& fhat_beta,
    int N_levels,
    int level) {

        // Get N from level in tree
        int M = T_alpha.rows();            // # of rows or columns in T
        int k = N_levels - level - 1;      // Reverse level ID
        int N = M / pow(2, k / 2 + 2);     // # of cells per leaf patch per side
        int B = N * pow(2, k / 2);         // # of rows or columns in block size

        // Create permutation matrix
        //    Alpha
        Vector<int> pi_alpha({0, 2, 3, 1});
        Matrix<double> P_alpha = blockPermutation<double>(pi_alpha, B);
        Matrix<double> P_alphaT = P_alpha.transpose();

        //    Beta
        Vector<int> pi_beta({1, 2, 3, 0});
        Matrix<double> P_beta = blockPermutation<double>(pi_beta, B);
        Matrix<double> P_betaT = P_beta.transpose();

        // Permutate DtN matricies
        Matrix<double> T_alpha_reordered = P_alpha * T_alpha;
        T_alpha_reordered = T_alpha_reordered * P_alphaT;
        Matrix<double> T_beta_reordered = P_beta * T_beta;
        T_beta_reordered = T_beta_reordered * P_betaT;

        // Permutate fhat vector
        Vector<double> fhat_alpha_reordered = P_alpha * fhat_alpha;
        Vector<double> fhat_beta_reordered = P_beta * fhat_beta;

        // Extract blocks
        Matrix<double> T_alpha_11 = T_alpha_reordered.extract(0*N, 0*N, 3*N, 3*N);
        Matrix<double> T_alpha_13 = T_alpha_reordered.extract(0*N, 3*N, 3*N, 1*N);
        Matrix<double> T_alpha_31 = T_alpha_reordered.extract(3*N, 0*N, 1*N, 3*N);
        Matrix<double> T_alpha_33 = T_alpha_reordered.extract(3*N, 3*N, 1*N, 1*N);

        Matrix<double> T_beta_22 = T_beta_reordered.extract(0*N, 0*N, 3*N, 3*N);
        Matrix<double> T_beta_23 = T_beta_reordered.extract(0*N, 3*N, 3*N, 1*N);
        Matrix<double> T_beta_32 = T_beta_reordered.extract(3*N, 0*N, 1*N, 3*N);
        Matrix<double> T_beta_33 = T_beta_reordered.extract(3*N, 3*N, 1*N, 1*N);

        Vector<double> fhat_alpha_1 = fhat_alpha_reordered.extract(0*N, 3*N);
        Vector<double> fhat_alpha_3 = fhat_alpha_reordered.extract(3*N, 1*N);
        Vector<double> fhat_beta_2 = fhat_beta_reordered.extract(0*N, 3*N);
        Vector<double> fhat_beta_3 = fhat_beta_reordered.extract(3*N, 1*N);

        // Perform merge operation
        Matrix<double> S = mergeOperationS(T_alpha_33, T_beta_33, T_alpha_31, T_beta_32);
        Matrix<double> T = mergeOperationT(T_alpha_11, T_beta_22, T_alpha_13, T_beta_23, S);
        Vector<double> fhat = mergeOperationF(T_alpha_13, T_beta_23, T_alpha_33, T_beta_33, fhat_alpha_1, fhat_beta_2, fhat_alpha_3, fhat_beta_3);
        Vector<double> w = mergeOperationW(T_alpha_33, T_beta_33, fhat_alpha_3, fhat_beta_3);

        Vector<int> pi_S({0, 3, 1, 4, 2, 5});
        Matrix<double> P_S = blockPermutation<double>(pi_S, B);
        Matrix<double> P_ST = P_S.transpose();
        S = S * P_ST;

        return merge_values{S, T, fhat, w};

}

/**
 * Performs a vertical merge.
 * @param  T_alpha    Alpha DtN matrix
 * @param  T_beta     Beta DtN matrix
 * @param  fhat_alpha Alpha fhat vector
 * @param  fhat_beta  Beta fhat vector
 * @param  N_levels   Number of levels in tree
 * @param  level      Current level in tree
 * @return            merge_values struct with parent patch data
 */
merge_values mergeVertical(
    Matrix<double>& T_alpha,
    Matrix<double>& T_beta,
    Vector<double>& fhat_alpha,
    Vector<double>& fhat_beta,
    int N_levels,
    int level) {

        // Get N from level in tree
        int M = T_alpha.rows();
        int k = N_levels - level - 1;
        int N = M / (pow(2, k/2 + 2) + pow(2, k/2 + 1));
        int B = N * pow(2, k / 2);         // # of rows or columns in block size

        // Create permutation matrix
        //    Alpha
        Vector<int> pi_alpha({0, 3, 1, 4, 2, 5});
        Matrix<double> P_alpha = blockPermutation<double>(pi_alpha, B);
        Matrix<double> P_alphaT = P_alpha.transpose();

        //    Beta
        Vector<int> pi_beta({0, 3, 2, 5, 1, 4});
        Matrix<double> P_beta = blockPermutation<double>(pi_beta, B);
        Matrix<double> P_betaT = P_beta.transpose();

        // Permutate DtN matricies
        Matrix<double> T_alpha_reordered = P_alpha * T_alpha;
        T_alpha_reordered = T_alpha_reordered * P_alphaT;
        Matrix<double> T_beta_reordered = P_beta * T_beta;
        T_beta_reordered = T_beta_reordered * P_betaT;

        // Permutate fhat vector
        Vector<double> fhat_alpha_reordered = P_alpha * fhat_alpha;
        Vector<double> fhat_beta_reordered = P_beta * fhat_beta;

        // Extract blocks
        Matrix<double> T_alpha_11 = T_alpha_reordered.extract(0*N, 0*N, 4*N, 4*N);
        Matrix<double> T_alpha_13 = T_alpha_reordered.extract(0*N, 4*N, 4*N, 2*N);
        Matrix<double> T_alpha_31 = T_alpha_reordered.extract(4*N, 0*N, 2*N, 4*N);
        Matrix<double> T_alpha_33 = T_alpha_reordered.extract(4*N, 4*N, 2*N, 2*N);

        Matrix<double> T_beta_22 = T_beta_reordered.extract(0*N, 0*N, 4*N, 4*N);
        Matrix<double> T_beta_23 = T_beta_reordered.extract(0*N, 4*N, 4*N, 2*N);
        Matrix<double> T_beta_32 = T_beta_reordered.extract(4*N, 0*N, 2*N, 4*N);
        Matrix<double> T_beta_33 = T_beta_reordered.extract(4*N, 4*N, 2*N, 2*N);

        Vector<double> fhat_alpha_1 = fhat_alpha_reordered.extract(0*N, 4*N);
        Vector<double> fhat_alpha_3 = fhat_alpha_reordered.extract(4*N, 2*N);
        Vector<double> fhat_beta_2 = fhat_beta_reordered.extract(0*N, 4*N);
        Vector<double> fhat_beta_3 = fhat_beta_reordered.extract(4*N, 2*N);

        // Perform merge operation
        Matrix<double> S = mergeOperationS(T_alpha_33, T_beta_33, T_alpha_31, T_beta_32);
        Matrix<double> T = mergeOperationT(T_alpha_11, T_beta_22, T_alpha_13, T_beta_23, S);
        Vector<double> fhat = mergeOperationF(T_alpha_13, T_beta_23, T_alpha_33, T_beta_33, fhat_alpha_1, fhat_beta_2, fhat_alpha_3, fhat_beta_3);
        Vector<double> w = mergeOperationW(T_alpha_33, T_beta_33, fhat_alpha_3, fhat_beta_3);

        // Reorder to return to WESN ordering
        Vector<int> pi_tau({0, 4, 1, 5, 2, 3, 6, 7});
        Matrix<double> P_tau = blockPermutation<double>(pi_tau, N);
        Matrix<double> P_tauT = P_tau.transpose();

        T = P_tau * T;
        T = T * P_tauT;
        fhat = P_tau * fhat;
        // w = P_tau * w;
        S = S * P_tauT;

        return merge_values{S, T, fhat, w};

}


} // END NAMESPACE hps
