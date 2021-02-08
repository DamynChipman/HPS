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
 * @TODO: Unit Test
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
 * Performs a horizontal merge for a provided DtN operator `T_level` at the specified `level`. Extracts the blocks corresponding to each side of the patch, intracts them into the merge operation matrices, and calls `mergeOperationS` and `mergeOperationT` to do the linear algebra.
 * @NOTE: This HAS been unit tested against Python and passed
 * @param T_level  Matrix for DtN operator at the specified level. Should be a (4*N x 4*N) matrix for each level away from the leaf level
 * @param N_levels Number of levels in the tree
 * @param level    Current level. Should be even for the horizontal merge
 * @return         Returns a std::pair of matrices (S,T)
 */
merge_values mergeHorizontal(const Matrix<double>& T_level, int N_levels, int level) {

    // Get number of entries in each section
    int N = (T_level.rows() / 4) * (N_levels - level);

    // Unpack input matrix into blocks
    Matrix<double> T_WW = T_level.extract(0*N, 0*N, N, N);
    Matrix<double> T_WE = T_level.extract(0*N, 1*N, N, N);
    Matrix<double> T_WS = T_level.extract(0*N, 2*N, N, N);
    Matrix<double> T_WN = T_level.extract(0*N, 3*N, N, N);

    Matrix<double> T_EW = T_level.extract(1*N, 0*N, N, N);
    Matrix<double> T_EE = T_level.extract(1*N, 1*N, N, N);
    Matrix<double> T_ES = T_level.extract(1*N, 2*N, N, N);
    Matrix<double> T_EN = T_level.extract(1*N, 3*N, N, N);

    Matrix<double> T_SW = T_level.extract(2*N, 0*N, N, N);
    Matrix<double> T_SE = T_level.extract(2*N, 1*N, N, N);
    Matrix<double> T_SS = T_level.extract(2*N, 2*N, N, N);
    Matrix<double> T_SN = T_level.extract(2*N, 3*N, N, N);

    Matrix<double> T_NW = T_level.extract(3*N, 0*N, N, N);
    Matrix<double> T_NE = T_level.extract(3*N, 1*N, N, N);
    Matrix<double> T_NS = T_level.extract(3*N, 2*N, N, N);
    Matrix<double> T_NN = T_level.extract(3*N, 3*N, N, N);

    // Create partitioned matrices
    Matrix<double> T_alpha_11(3*N, 3*N);
    T_alpha_11.intract(0*N, 0*N, T_WW);
    T_alpha_11.intract(0*N, 1*N, T_WS);
    T_alpha_11.intract(0*N, 2*N, T_WN);
    T_alpha_11.intract(1*N, 0*N, T_SW);
    T_alpha_11.intract(1*N, 1*N, T_SS);
    T_alpha_11.intract(1*N, 2*N, T_SN);
    T_alpha_11.intract(2*N, 0*N, T_NW);
    T_alpha_11.intract(2*N, 1*N, T_NS);
    T_alpha_11.intract(2*N, 2*N, T_NN);

    Matrix<double> T_alpha_13(3*N, 1*N);
    T_alpha_13.intract(0*N, 0*N, T_WE);
    T_alpha_13.intract(1*N, 0*N, T_SE);
    T_alpha_13.intract(2*N, 0*N, T_NE);

    Matrix<double> T_alpha_31(1*N, 3*N);
    T_alpha_31.intract(0*N, 0*N, T_EW);
    T_alpha_31.intract(0*N, 1*N, T_ES);
    T_alpha_31.intract(0*N, 2*N, T_EN);

    Matrix<double> T_alpha_33(1*N, 1*N);
    T_alpha_33.intract(0*N, 0*N, T_EE);

    Matrix<double> T_beta_22(3*N, 3*N);
    T_beta_22.intract(0*N, 0*N, T_EE);
    T_beta_22.intract(0*N, 1*N, T_ES);
    T_beta_22.intract(0*N, 2*N, T_EN);
    T_beta_22.intract(1*N, 0*N, T_SE);
    T_beta_22.intract(1*N, 1*N, T_SS);
    T_beta_22.intract(1*N, 2*N, T_SN);
    T_beta_22.intract(2*N, 0*N, T_NE);
    T_beta_22.intract(2*N, 1*N, T_NS);
    T_beta_22.intract(2*N, 2*N, T_NN);

    Matrix<double> T_beta_23(3*N, 1*N);
    T_beta_23.intract(0*N, 0*N, T_EW);
    T_beta_23.intract(1*N, 0*N, T_SW);
    T_beta_23.intract(2*N, 0*N, T_NW);

    Matrix<double> T_beta_32(1*N, 3*N);
    T_beta_32.intract(0*N, 0*N, T_WE);
    T_beta_32.intract(0*N, 1*N, T_WS);
    T_beta_32.intract(0*N, 2*N, T_WN);

    Matrix<double> T_beta_33(1*N, 1*N);
    T_beta_33.intract(0*N, 0*N, T_WW);

    // Perform merge operation
    Matrix<double> S = mergeOperationS(T_alpha_33, T_beta_33, T_alpha_31, T_beta_32);
    Matrix<double> T = mergeOperationT(T_alpha_11, T_beta_22, T_alpha_13, T_beta_23, S);

    return merge_values{S, T};
}

/**
 * Performs a vertical merge for a provided DtN operator `T_level` at the specified `level`. Extracts the blocks corresponding to each side of the patch, intracts them into the merge operation matrices, and calls `mergeOperationS` and `mergeOperationT` to do the linear algebra. After the linear algebra, reorders S and T back into the WESN ordering for a square (via a permutation matrix multiplication).
 * @NOTE: This HAS been unit tested against Python and passed
 * @param T_level  Matrix for DtN operator at the specified level. Should be a (6*N x 6*N) matrix for each level away from the leaf level
 * @param N_levels Number of levels in the tree
 * @param level    Current level. Should be odd for the vertical merge
 * @return         Returns a std::pair of matrices (S,T)
 */
merge_values mergeVertical(const Matrix<double>& T_level, int N_levels, int level) {

    // Get number of entries in each section
    int N = (T_level.rows() / 6) * (N_levels - level - 1);

    // Unpack input matrix into blocks
    Matrix<double> T_WW_aa = T_level.extract(0*N, 0*N, N, N);
    Matrix<double> T_WS_aa = T_level.extract(0*N, 1*N, N, N);
    Matrix<double> T_WN_aa = T_level.extract(0*N, 2*N, N, N);
    Matrix<double> T_WE_ab = T_level.extract(0*N, 3*N, N, N);
    Matrix<double> T_WS_ab = T_level.extract(0*N, 4*N, N, N);
    Matrix<double> T_WN_ab = T_level.extract(0*N, 5*N, N, N);

    Matrix<double> T_SW_aa = T_level.extract(1*N, 0*N, N, N);
    Matrix<double> T_SS_aa = T_level.extract(1*N, 1*N, N, N);
    Matrix<double> T_SN_aa = T_level.extract(1*N, 2*N, N, N);
    Matrix<double> T_SE_ab = T_level.extract(1*N, 3*N, N, N);
    Matrix<double> T_SS_ab = T_level.extract(1*N, 4*N, N, N);
    Matrix<double> T_SN_ab = T_level.extract(1*N, 5*N, N, N);

    Matrix<double> T_NW_aa = T_level.extract(2*N, 0*N, N, N);
    Matrix<double> T_NS_aa = T_level.extract(2*N, 1*N, N, N);
    Matrix<double> T_NN_aa = T_level.extract(2*N, 2*N, N, N);
    Matrix<double> T_NE_ab = T_level.extract(2*N, 3*N, N, N);
    Matrix<double> T_NS_ab = T_level.extract(2*N, 4*N, N, N);
    Matrix<double> T_NN_ab = T_level.extract(2*N, 5*N, N, N);

    Matrix<double> T_EW_ba = T_level.extract(3*N, 0*N, N, N);
    Matrix<double> T_ES_ba = T_level.extract(3*N, 1*N, N, N);
    Matrix<double> T_EN_ba = T_level.extract(3*N, 2*N, N, N);
    Matrix<double> T_EE_bb = T_level.extract(3*N, 3*N, N, N);
    Matrix<double> T_ES_bb = T_level.extract(3*N, 4*N, N, N);
    Matrix<double> T_EN_bb = T_level.extract(3*N, 5*N, N, N);

    Matrix<double> T_SW_ba = T_level.extract(4*N, 0*N, N, N);
    Matrix<double> T_SS_ba = T_level.extract(4*N, 1*N, N, N);
    Matrix<double> T_SN_ba = T_level.extract(4*N, 2*N, N, N);
    Matrix<double> T_SE_bb = T_level.extract(4*N, 3*N, N, N);
    Matrix<double> T_SS_bb = T_level.extract(4*N, 4*N, N, N);
    Matrix<double> T_SN_bb = T_level.extract(4*N, 5*N, N, N);

    Matrix<double> T_NW_ba = T_level.extract(5*N, 0*N, N, N);
    Matrix<double> T_NS_ba = T_level.extract(5*N, 1*N, N, N);
    Matrix<double> T_NN_ba = T_level.extract(5*N, 2*N, N, N);
    Matrix<double> T_NE_bb = T_level.extract(5*N, 3*N, N, N);
    Matrix<double> T_NS_bb = T_level.extract(5*N, 4*N, N, N);
    Matrix<double> T_NN_bb = T_level.extract(5*N, 5*N, N, N);

    // Create partitioned matricies
    Matrix<double> T_alpha_11(4*N, 4*N);
    T_alpha_11.intract(0*N, 0*N, T_WW_aa);
    T_alpha_11.intract(0*N, 1*N, T_WE_ab);
    T_alpha_11.intract(0*N, 2*N, T_WS_aa);
    T_alpha_11.intract(0*N, 3*N, T_WS_ab);
    T_alpha_11.intract(1*N, 0*N, T_EW_ba);
    T_alpha_11.intract(1*N, 1*N, T_EE_bb);
    T_alpha_11.intract(1*N, 2*N, T_ES_ba);
    T_alpha_11.intract(1*N, 3*N, T_ES_bb);
    T_alpha_11.intract(2*N, 0*N, T_SW_aa);
    T_alpha_11.intract(2*N, 1*N, T_SE_ab);
    T_alpha_11.intract(2*N, 2*N, T_SS_aa);
    T_alpha_11.intract(2*N, 3*N, T_SS_ab);
    T_alpha_11.intract(3*N, 0*N, T_SW_ba);
    T_alpha_11.intract(3*N, 1*N, T_SE_bb);
    T_alpha_11.intract(3*N, 2*N, T_SS_ba);
    T_alpha_11.intract(3*N, 3*N, T_SS_bb);

    Matrix<double> T_alpha_13(4*N, 2*N);
    T_alpha_13.intract(0*N, 0*N, T_WN_aa);
    T_alpha_13.intract(0*N, 1*N, T_WN_ab);
    T_alpha_13.intract(1*N, 0*N, T_EN_ba);
    T_alpha_13.intract(1*N, 1*N, T_EN_bb);
    T_alpha_13.intract(2*N, 0*N, T_SN_aa);
    T_alpha_13.intract(2*N, 1*N, T_SN_ab);
    T_alpha_13.intract(3*N, 0*N, T_SN_ba);
    T_alpha_13.intract(3*N, 1*N, T_SN_bb);

    Matrix<double> T_alpha_31(2*N, 4*N);
    T_alpha_31.intract(0*N, 0*N, T_NW_aa);
    T_alpha_31.intract(0*N, 1*N, T_NE_ab);
    T_alpha_31.intract(0*N, 2*N, T_NS_aa);
    T_alpha_31.intract(0*N, 3*N, T_NS_ab);
    T_alpha_31.intract(1*N, 0*N, T_NW_ba);
    T_alpha_31.intract(1*N, 1*N, T_NE_bb);
    T_alpha_31.intract(1*N, 2*N, T_NS_ba);
    T_alpha_31.intract(1*N, 3*N, T_NS_bb);

    Matrix<double> T_alpha_33(2*N, 2*N);
    T_alpha_33.intract(0*N, 0*N, T_NN_aa);
    T_alpha_33.intract(0*N, 1*N, T_NN_ab);
    T_alpha_33.intract(1*N, 0*N, T_NN_ba);
    T_alpha_33.intract(1*N, 1*N, T_NN_bb);

    Matrix<double> T_beta_22(4*N, 4*N);
    T_beta_22.intract(0*N, 0*N, T_WW_aa);
    T_beta_22.intract(0*N, 1*N, T_WE_ab);
    T_beta_22.intract(0*N, 2*N, T_WN_aa);
    T_beta_22.intract(0*N, 3*N, T_WN_ab);
    T_beta_22.intract(1*N, 0*N, T_EW_ba);
    T_beta_22.intract(1*N, 1*N, T_EE_bb);
    T_beta_22.intract(1*N, 2*N, T_EN_ba);
    T_beta_22.intract(1*N, 3*N, T_EN_bb);
    T_beta_22.intract(2*N, 0*N, T_NW_aa);
    T_beta_22.intract(2*N, 1*N, T_NE_ab);
    T_beta_22.intract(2*N, 2*N, T_NN_aa);
    T_beta_22.intract(2*N, 3*N, T_NN_ab);
    T_beta_22.intract(3*N, 0*N, T_NW_ba);
    T_beta_22.intract(3*N, 1*N, T_NE_bb);
    T_beta_22.intract(3*N, 2*N, T_NN_ba);
    T_beta_22.intract(3*N, 3*N, T_NN_bb);

    Matrix<double> T_beta_23(4*N, 2*N);
    T_beta_23.intract(0*N, 0*N, T_WS_aa);
    T_beta_23.intract(0*N, 1*N, T_WS_ab);
    T_beta_23.intract(1*N, 0*N, T_ES_ba);
    T_beta_23.intract(1*N, 1*N, T_ES_bb);
    T_beta_23.intract(2*N, 0*N, T_NS_aa);
    T_beta_23.intract(2*N, 1*N, T_NS_ab);
    T_beta_23.intract(3*N, 0*N, T_NS_ba);
    T_beta_23.intract(3*N, 1*N, T_NS_bb);

    Matrix<double> T_beta_32(2*N, 4*N);
    T_beta_32.intract(0*N, 0*N, T_SW_aa);
    T_beta_32.intract(0*N, 1*N, T_SE_ab);
    T_beta_32.intract(0*N, 2*N, T_SN_aa);
    T_beta_32.intract(0*N, 3*N, T_SN_ab);
    T_beta_32.intract(1*N, 0*N, T_SW_ba);
    T_beta_32.intract(1*N, 1*N, T_SE_bb);
    T_beta_32.intract(1*N, 2*N, T_SN_ba);
    T_beta_32.intract(1*N, 3*N, T_SN_bb);

    Matrix<double> T_beta_33(2*N, 2*N);
    T_beta_33.intract(0*N, 0*N, T_SS_aa);
    T_beta_33.intract(0*N, 1*N, T_SS_ab);
    T_beta_33.intract(1*N, 0*N, T_SS_ba);
    T_beta_33.intract(1*N, 1*N, T_SS_bb);

    // Perform merge operation
    Matrix<double> S = mergeOperationS(T_alpha_33, T_beta_33, T_alpha_31, T_beta_32);
    Matrix<double> T = mergeOperationT(T_alpha_11, T_beta_22, T_alpha_13, T_beta_23, S);

    // Reorder for square WESN ordering
    Matrix<double> reorder_matrix(8*N, 8*N);
    double z = 0.0;
    Matrix<double> eyeN(N, N, z);
    eyeN.identify();
    reorder_matrix.intract(0*N, 4*N, eyeN);
    reorder_matrix.intract(1*N, 0*N, eyeN);
    reorder_matrix.intract(2*N, 5*N, eyeN);
    reorder_matrix.intract(3*N, 1*N, eyeN);
    reorder_matrix.intract(4*N, 2*N, eyeN);
    reorder_matrix.intract(5*N, 3*N, eyeN);
    reorder_matrix.intract(6*N, 6*N, eyeN);
    reorder_matrix.intract(7*N, 7*N, eyeN);
    T = reorder_matrix * T;

    // @@TODO: Rerorder S matrix (if necessary)


    return merge_values{S, T};
}


merge_values mergeHorizontal(
    const Matrix<double>& T_alpha,
    const Matrix<double>& T_beta,
    const Vector<double>& fhat_alpha,
    const Vector<double>& fhat_beta,
    int N_levels,
    int level) {

        // Get number of entries in each section
        int N = (T_alpha.rows() / 4) * (N_levels - level);

        // Unpack input data into blocks
        //   Alpha Patch
        Matrix<double> T_WW_alpha = T_alpha.extract(0*N, 0*N, N, N);
        Matrix<double> T_WE_alpha = T_alpha.extract(0*N, 1*N, N, N);
        Matrix<double> T_WS_alpha = T_alpha.extract(0*N, 2*N, N, N);
        Matrix<double> T_WN_alpha = T_alpha.extract(0*N, 3*N, N, N);

        Matrix<double> T_EW_alpha = T_alpha.extract(1*N, 0*N, N, N);
        Matrix<double> T_EE_alpha = T_alpha.extract(1*N, 1*N, N, N);
        Matrix<double> T_ES_alpha = T_alpha.extract(1*N, 2*N, N, N);
        Matrix<double> T_EN_alpha = T_alpha.extract(1*N, 3*N, N, N);

        Matrix<double> T_SW_alpha = T_alpha.extract(2*N, 0*N, N, N);
        Matrix<double> T_SE_alpha = T_alpha.extract(2*N, 1*N, N, N);
        Matrix<double> T_SS_alpha = T_alpha.extract(2*N, 2*N, N, N);
        Matrix<double> T_SN_alpha = T_alpha.extract(2*N, 3*N, N, N);

        Matrix<double> T_NW_alpha = T_alpha.extract(3*N, 0*N, N, N);
        Matrix<double> T_NE_alpha = T_alpha.extract(3*N, 1*N, N, N);
        Matrix<double> T_NS_alpha = T_alpha.extract(3*N, 2*N, N, N);
        Matrix<double> T_NN_alpha = T_alpha.extract(3*N, 3*N, N, N);

        Vector<double> fhat_W_alpha = fhat_alpha.extract(0*N, N);
        Vector<double> fhat_E_alpha = fhat_alpha.extract(1*N, N);
        Vector<double> fhat_S_alpha = fhat_alpha.extract(2*N, N);
        Vector<double> fhat_N_alpha = fhat_alpha.extract(3*N, N);

        //   Beta Patch
        Matrix<double> T_WW_beta = T_beta.extract(0*N, 0*N, N, N);
        Matrix<double> T_WE_beta = T_beta.extract(0*N, 1*N, N, N);
        Matrix<double> T_WS_beta = T_beta.extract(0*N, 2*N, N, N);
        Matrix<double> T_WN_beta = T_beta.extract(0*N, 3*N, N, N);

        Matrix<double> T_EW_beta = T_beta.extract(1*N, 0*N, N, N);
        Matrix<double> T_EE_beta = T_beta.extract(1*N, 1*N, N, N);
        Matrix<double> T_ES_beta = T_beta.extract(1*N, 2*N, N, N);
        Matrix<double> T_EN_beta = T_beta.extract(1*N, 3*N, N, N);

        Matrix<double> T_SW_beta = T_beta.extract(2*N, 0*N, N, N);
        Matrix<double> T_SE_beta = T_beta.extract(2*N, 1*N, N, N);
        Matrix<double> T_SS_beta = T_beta.extract(2*N, 2*N, N, N);
        Matrix<double> T_SN_beta = T_beta.extract(2*N, 3*N, N, N);

        Matrix<double> T_NW_beta = T_beta.extract(3*N, 0*N, N, N);
        Matrix<double> T_NE_beta = T_beta.extract(3*N, 1*N, N, N);
        Matrix<double> T_NS_beta = T_beta.extract(3*N, 2*N, N, N);
        Matrix<double> T_NN_beta = T_beta.extract(3*N, 3*N, N, N);

        Vector<double> fhat_W_beta = fhat_beta.extract(0*N, N);
        Vector<double> fhat_E_beta = fhat_beta.extract(1*N, N);
        Vector<double> fhat_S_beta = fhat_beta.extract(2*N, N);
        Vector<double> fhat_N_beta = fhat_beta.extract(3*N, N);

        // Create partitioned matrices
        Matrix<double> T_alpha_11(3*N, 3*N);
        T_alpha_11.intract(0*N, 0*N, T_WW_alpha);
        T_alpha_11.intract(0*N, 1*N, T_WS_alpha);
        T_alpha_11.intract(0*N, 2*N, T_WN_alpha);
        T_alpha_11.intract(1*N, 0*N, T_SW_alpha);
        T_alpha_11.intract(1*N, 1*N, T_SS_alpha);
        T_alpha_11.intract(1*N, 2*N, T_SN_alpha);
        T_alpha_11.intract(2*N, 0*N, T_NW_alpha);
        T_alpha_11.intract(2*N, 1*N, T_NS_alpha);
        T_alpha_11.intract(2*N, 2*N, T_NN_alpha);

        Matrix<double> T_alpha_13(3*N, 1*N);
        T_alpha_13.intract(0*N, 0*N, T_WE_alpha);
        T_alpha_13.intract(1*N, 0*N, T_SE_alpha);
        T_alpha_13.intract(2*N, 0*N, T_NE_alpha);

        Matrix<double> T_alpha_31(1*N, 3*N);
        T_alpha_31.intract(0*N, 0*N, T_EW_alpha);
        T_alpha_31.intract(0*N, 1*N, T_ES_alpha);
        T_alpha_31.intract(0*N, 2*N, T_EN_alpha);

        Matrix<double> T_alpha_33(1*N, 1*N);
        T_alpha_33.intract(0*N, 0*N, T_EE_alpha);

        Matrix<double> T_beta_22(3*N, 3*N);
        T_beta_22.intract(0*N, 0*N, T_EE_beta);
        T_beta_22.intract(0*N, 1*N, T_ES_beta);
        T_beta_22.intract(0*N, 2*N, T_EN_beta);
        T_beta_22.intract(1*N, 0*N, T_SE_beta);
        T_beta_22.intract(1*N, 1*N, T_SS_beta);
        T_beta_22.intract(1*N, 2*N, T_SN_beta);
        T_beta_22.intract(2*N, 0*N, T_NE_beta);
        T_beta_22.intract(2*N, 1*N, T_NS_beta);
        T_beta_22.intract(2*N, 2*N, T_NN_beta);

        Matrix<double> T_beta_23(3*N, 1*N);
        T_beta_23.intract(0*N, 0*N, T_EW_beta);
        T_beta_23.intract(1*N, 0*N, T_SW_beta);
        T_beta_23.intract(2*N, 0*N, T_NW_beta);

        Matrix<double> T_beta_32(1*N, 3*N);
        T_beta_32.intract(0*N, 0*N, T_WE_beta);
        T_beta_32.intract(0*N, 1*N, T_WS_beta);
        T_beta_32.intract(0*N, 2*N, T_WN_beta);

        Matrix<double> T_beta_33(1*N, 1*N);
        T_beta_33.intract(0*N, 0*N, T_WW_beta);

        Vector<double> fhat_alpha_1(3*N);
        fhat_alpha_1.intract(0*N, fhat_W_alpha);
        fhat_alpha_1.intract(1*N, fhat_S_alpha);
        fhat_alpha_1.intract(2*N, fhat_N_alpha);

        Vector<double> fhat_alpha_3(1*N);
        fhat_alpha_3.intract(0*N, fhat_E_alpha);

        Vector<double> fhat_beta_2(3*N);
        fhat_beta_2.intract(0*N, fhat_E_beta);
        fhat_beta_2.intract(1*N, fhat_S_beta);
        fhat_beta_2.intract(2*N, fhat_N_beta);

        Vector<double> fhat_beta_3(1*N);
        fhat_beta_3.intract(0*N, fhat_W_beta);

        // Perform merge operation
        Matrix<double> S = mergeOperationS(T_alpha_33, T_beta_33, T_alpha_31, T_beta_32);
        Matrix<double> T = mergeOperationT(T_alpha_11, T_beta_22, T_alpha_13, T_beta_23, S);
        Vector<double> fhat = mergeOperationF(T_alpha_13, T_beta_23, T_alpha_33, T_beta_33, fhat_alpha_1, fhat_beta_2, fhat_alpha_3, fhat_beta_3);

        return merge_values{S, T, fhat};

}


merge_values mergeVertical(
    const Matrix<double>& T_alpha,
    const Matrix<double>& T_beta,
    const Vector<double>& fhat_alpha,
    const Vector<double>& fhat_beta,
    int N_levels,
    int level) {

        // Get number of entries in each section
        int N = (T_alpha.rows() / 6) * (N_levels - level - 1);

        // Unpack input matrix into blocks
        //   Alpha Patch
        Matrix<double> T_WW_aa_alpha = T_alpha.extract(0*N, 0*N, N, N);
        Matrix<double> T_WS_aa_alpha = T_alpha.extract(0*N, 1*N, N, N);
        Matrix<double> T_WN_aa_alpha = T_alpha.extract(0*N, 2*N, N, N);
        Matrix<double> T_WE_ab_alpha = T_alpha.extract(0*N, 3*N, N, N);
        Matrix<double> T_WS_ab_alpha = T_alpha.extract(0*N, 4*N, N, N);
        Matrix<double> T_WN_ab_alpha = T_alpha.extract(0*N, 5*N, N, N);

        Matrix<double> T_SW_aa_alpha = T_alpha.extract(1*N, 0*N, N, N);
        Matrix<double> T_SS_aa_alpha = T_alpha.extract(1*N, 1*N, N, N);
        Matrix<double> T_SN_aa_alpha = T_alpha.extract(1*N, 2*N, N, N);
        Matrix<double> T_SE_ab_alpha = T_alpha.extract(1*N, 3*N, N, N);
        Matrix<double> T_SS_ab_alpha = T_alpha.extract(1*N, 4*N, N, N);
        Matrix<double> T_SN_ab_alpha = T_alpha.extract(1*N, 5*N, N, N);

        Matrix<double> T_NW_aa_alpha = T_alpha.extract(2*N, 0*N, N, N);
        Matrix<double> T_NS_aa_alpha = T_alpha.extract(2*N, 1*N, N, N);
        Matrix<double> T_NN_aa_alpha = T_alpha.extract(2*N, 2*N, N, N);
        Matrix<double> T_NE_ab_alpha = T_alpha.extract(2*N, 3*N, N, N);
        Matrix<double> T_NS_ab_alpha = T_alpha.extract(2*N, 4*N, N, N);
        Matrix<double> T_NN_ab_alpha = T_alpha.extract(2*N, 5*N, N, N);

        Matrix<double> T_EW_ba_alpha = T_alpha.extract(3*N, 0*N, N, N);
        Matrix<double> T_ES_ba_alpha = T_alpha.extract(3*N, 1*N, N, N);
        Matrix<double> T_EN_ba_alpha = T_alpha.extract(3*N, 2*N, N, N);
        Matrix<double> T_EE_bb_alpha = T_alpha.extract(3*N, 3*N, N, N);
        Matrix<double> T_ES_bb_alpha = T_alpha.extract(3*N, 4*N, N, N);
        Matrix<double> T_EN_bb_alpha = T_alpha.extract(3*N, 5*N, N, N);

        Matrix<double> T_SW_ba_alpha = T_alpha.extract(4*N, 0*N, N, N);
        Matrix<double> T_SS_ba_alpha = T_alpha.extract(4*N, 1*N, N, N);
        Matrix<double> T_SN_ba_alpha = T_alpha.extract(4*N, 2*N, N, N);
        Matrix<double> T_SE_bb_alpha = T_alpha.extract(4*N, 3*N, N, N);
        Matrix<double> T_SS_bb_alpha = T_alpha.extract(4*N, 4*N, N, N);
        Matrix<double> T_SN_bb_alpha = T_alpha.extract(4*N, 5*N, N, N);

        Matrix<double> T_NW_ba_alpha = T_alpha.extract(5*N, 0*N, N, N);
        Matrix<double> T_NS_ba_alpha = T_alpha.extract(5*N, 1*N, N, N);
        Matrix<double> T_NN_ba_alpha = T_alpha.extract(5*N, 2*N, N, N);
        Matrix<double> T_NE_bb_alpha = T_alpha.extract(5*N, 3*N, N, N);
        Matrix<double> T_NS_bb_alpha = T_alpha.extract(5*N, 4*N, N, N);
        Matrix<double> T_NN_bb_alpha = T_alpha.extract(5*N, 5*N, N, N);

        Vector<double> fhat_alpha_aW = fhat_alpha.extract(0*N, N);
        Vector<double> fhat_alpha_aS = fhat_alpha.extract(1*N, N);
        Vector<double> fhat_alpha_aN = fhat_alpha.extract(2*N, N);
        Vector<double> fhat_alpha_bE = fhat_alpha.extract(3*N, N);
        Vector<double> fhat_alpha_bS = fhat_alpha.extract(4*N, N);
        Vector<double> fhat_alpha_bN = fhat_alpha.extract(5*N, N);

        //   Beta Patch
        Matrix<double> T_WW_aa_beta = T_beta.extract(0*N, 0*N, N, N);
        Matrix<double> T_WS_aa_beta = T_beta.extract(0*N, 1*N, N, N);
        Matrix<double> T_WN_aa_beta = T_beta.extract(0*N, 2*N, N, N);
        Matrix<double> T_WE_ab_beta = T_beta.extract(0*N, 3*N, N, N);
        Matrix<double> T_WS_ab_beta = T_beta.extract(0*N, 4*N, N, N);
        Matrix<double> T_WN_ab_beta = T_beta.extract(0*N, 5*N, N, N);

        Matrix<double> T_SW_aa_beta = T_beta.extract(1*N, 0*N, N, N);
        Matrix<double> T_SS_aa_beta = T_beta.extract(1*N, 1*N, N, N);
        Matrix<double> T_SN_aa_beta = T_beta.extract(1*N, 2*N, N, N);
        Matrix<double> T_SE_ab_beta = T_beta.extract(1*N, 3*N, N, N);
        Matrix<double> T_SS_ab_beta = T_beta.extract(1*N, 4*N, N, N);
        Matrix<double> T_SN_ab_beta = T_beta.extract(1*N, 5*N, N, N);

        Matrix<double> T_NW_aa_beta = T_beta.extract(2*N, 0*N, N, N);
        Matrix<double> T_NS_aa_beta = T_beta.extract(2*N, 1*N, N, N);
        Matrix<double> T_NN_aa_beta = T_beta.extract(2*N, 2*N, N, N);
        Matrix<double> T_NE_ab_beta = T_beta.extract(2*N, 3*N, N, N);
        Matrix<double> T_NS_ab_beta = T_beta.extract(2*N, 4*N, N, N);
        Matrix<double> T_NN_ab_beta = T_beta.extract(2*N, 5*N, N, N);

        Matrix<double> T_EW_ba_beta = T_beta.extract(3*N, 0*N, N, N);
        Matrix<double> T_ES_ba_beta = T_beta.extract(3*N, 1*N, N, N);
        Matrix<double> T_EN_ba_beta = T_beta.extract(3*N, 2*N, N, N);
        Matrix<double> T_EE_bb_beta = T_beta.extract(3*N, 3*N, N, N);
        Matrix<double> T_ES_bb_beta = T_beta.extract(3*N, 4*N, N, N);
        Matrix<double> T_EN_bb_beta = T_beta.extract(3*N, 5*N, N, N);

        Matrix<double> T_SW_ba_beta = T_beta.extract(4*N, 0*N, N, N);
        Matrix<double> T_SS_ba_beta = T_beta.extract(4*N, 1*N, N, N);
        Matrix<double> T_SN_ba_beta = T_beta.extract(4*N, 2*N, N, N);
        Matrix<double> T_SE_bb_beta = T_beta.extract(4*N, 3*N, N, N);
        Matrix<double> T_SS_bb_beta = T_beta.extract(4*N, 4*N, N, N);
        Matrix<double> T_SN_bb_beta = T_beta.extract(4*N, 5*N, N, N);

        Matrix<double> T_NW_ba_beta = T_beta.extract(5*N, 0*N, N, N);
        Matrix<double> T_NS_ba_beta = T_beta.extract(5*N, 1*N, N, N);
        Matrix<double> T_NN_ba_beta = T_beta.extract(5*N, 2*N, N, N);
        Matrix<double> T_NE_bb_beta = T_beta.extract(5*N, 3*N, N, N);
        Matrix<double> T_NS_bb_beta = T_beta.extract(5*N, 4*N, N, N);
        Matrix<double> T_NN_bb_beta = T_beta.extract(5*N, 5*N, N, N);

        Vector<double> fhat_beta_aW = fhat_beta.extract(0*N, N);
        Vector<double> fhat_beta_aS = fhat_beta.extract(1*N, N);
        Vector<double> fhat_beta_aN = fhat_beta.extract(2*N, N);
        Vector<double> fhat_beta_bE = fhat_beta.extract(3*N, N);
        Vector<double> fhat_beta_bS = fhat_beta.extract(4*N, N);
        Vector<double> fhat_beta_bN = fhat_beta.extract(5*N, N);

        // Create partitioned matricies
        Matrix<double> T_alpha_11(4*N, 4*N);
        T_alpha_11.intract(0*N, 0*N, T_WW_aa_alpha);
        T_alpha_11.intract(0*N, 1*N, T_WE_ab_alpha);
        T_alpha_11.intract(0*N, 2*N, T_WS_aa_alpha);
        T_alpha_11.intract(0*N, 3*N, T_WS_ab_alpha);
        T_alpha_11.intract(1*N, 0*N, T_EW_ba_alpha);
        T_alpha_11.intract(1*N, 1*N, T_EE_bb_alpha);
        T_alpha_11.intract(1*N, 2*N, T_ES_ba_alpha);
        T_alpha_11.intract(1*N, 3*N, T_ES_bb_alpha);
        T_alpha_11.intract(2*N, 0*N, T_SW_aa_alpha);
        T_alpha_11.intract(2*N, 1*N, T_SE_ab_alpha);
        T_alpha_11.intract(2*N, 2*N, T_SS_aa_alpha);
        T_alpha_11.intract(2*N, 3*N, T_SS_ab_alpha);
        T_alpha_11.intract(3*N, 0*N, T_SW_ba_alpha);
        T_alpha_11.intract(3*N, 1*N, T_SE_bb_alpha);
        T_alpha_11.intract(3*N, 2*N, T_SS_ba_alpha);
        T_alpha_11.intract(3*N, 3*N, T_SS_bb_alpha);

        Matrix<double> T_alpha_13(4*N, 2*N);
        T_alpha_13.intract(0*N, 0*N, T_WN_aa_alpha);
        T_alpha_13.intract(0*N, 1*N, T_WN_ab_alpha);
        T_alpha_13.intract(1*N, 0*N, T_EN_ba_alpha);
        T_alpha_13.intract(1*N, 1*N, T_EN_bb_alpha);
        T_alpha_13.intract(2*N, 0*N, T_SN_aa_alpha);
        T_alpha_13.intract(2*N, 1*N, T_SN_ab_alpha);
        T_alpha_13.intract(3*N, 0*N, T_SN_ba_alpha);
        T_alpha_13.intract(3*N, 1*N, T_SN_bb_alpha);

        Matrix<double> T_alpha_31(2*N, 4*N);
        T_alpha_31.intract(0*N, 0*N, T_NW_aa_alpha);
        T_alpha_31.intract(0*N, 1*N, T_NE_ab_alpha);
        T_alpha_31.intract(0*N, 2*N, T_NS_aa_alpha);
        T_alpha_31.intract(0*N, 3*N, T_NS_ab_alpha);
        T_alpha_31.intract(1*N, 0*N, T_NW_ba_alpha);
        T_alpha_31.intract(1*N, 1*N, T_NE_bb_alpha);
        T_alpha_31.intract(1*N, 2*N, T_NS_ba_alpha);
        T_alpha_31.intract(1*N, 3*N, T_NS_bb_alpha);

        Matrix<double> T_alpha_33(2*N, 2*N);
        T_alpha_33.intract(0*N, 0*N, T_NN_aa_alpha);
        T_alpha_33.intract(0*N, 1*N, T_NN_ab_alpha);
        T_alpha_33.intract(1*N, 0*N, T_NN_ba_alpha);
        T_alpha_33.intract(1*N, 1*N, T_NN_bb_alpha);

        Matrix<double> T_beta_22(4*N, 4*N);
        T_beta_22.intract(0*N, 0*N, T_WW_aa_beta);
        T_beta_22.intract(0*N, 1*N, T_WE_ab_beta);
        T_beta_22.intract(0*N, 2*N, T_WN_aa_beta);
        T_beta_22.intract(0*N, 3*N, T_WN_ab_beta);
        T_beta_22.intract(1*N, 0*N, T_EW_ba_beta);
        T_beta_22.intract(1*N, 1*N, T_EE_bb_beta);
        T_beta_22.intract(1*N, 2*N, T_EN_ba_beta);
        T_beta_22.intract(1*N, 3*N, T_EN_bb_beta);
        T_beta_22.intract(2*N, 0*N, T_NW_aa_beta);
        T_beta_22.intract(2*N, 1*N, T_NE_ab_beta);
        T_beta_22.intract(2*N, 2*N, T_NN_aa_beta);
        T_beta_22.intract(2*N, 3*N, T_NN_ab_beta);
        T_beta_22.intract(3*N, 0*N, T_NW_ba_beta);
        T_beta_22.intract(3*N, 1*N, T_NE_bb_beta);
        T_beta_22.intract(3*N, 2*N, T_NN_ba_beta);
        T_beta_22.intract(3*N, 3*N, T_NN_bb_beta);

        Matrix<double> T_beta_23(4*N, 2*N);
        T_beta_23.intract(0*N, 0*N, T_WS_aa_beta);
        T_beta_23.intract(0*N, 1*N, T_WS_ab_beta);
        T_beta_23.intract(1*N, 0*N, T_ES_ba_beta);
        T_beta_23.intract(1*N, 1*N, T_ES_bb_beta);
        T_beta_23.intract(2*N, 0*N, T_NS_aa_beta);
        T_beta_23.intract(2*N, 1*N, T_NS_ab_beta);
        T_beta_23.intract(3*N, 0*N, T_NS_ba_beta);
        T_beta_23.intract(3*N, 1*N, T_NS_bb_beta);

        Matrix<double> T_beta_32(2*N, 4*N);
        T_beta_32.intract(0*N, 0*N, T_SW_aa_beta);
        T_beta_32.intract(0*N, 1*N, T_SE_ab_beta);
        T_beta_32.intract(0*N, 2*N, T_SN_aa_beta);
        T_beta_32.intract(0*N, 3*N, T_SN_ab_beta);
        T_beta_32.intract(1*N, 0*N, T_SW_ba_beta);
        T_beta_32.intract(1*N, 1*N, T_SE_bb_beta);
        T_beta_32.intract(1*N, 2*N, T_SN_ba_beta);
        T_beta_32.intract(1*N, 3*N, T_SN_bb_beta);

        Matrix<double> T_beta_33(2*N, 2*N);
        T_beta_33.intract(0*N, 0*N, T_SS_aa_beta);
        T_beta_33.intract(0*N, 1*N, T_SS_ab_beta);
        T_beta_33.intract(1*N, 0*N, T_SS_ba_beta);
        T_beta_33.intract(1*N, 1*N, T_SS_bb_beta);

        Vector<double> fhat_alpha_1(4*N);
        fhat_alpha_1.intract(0*N, fhat_alpha_aW);
        fhat_alpha_1.intract(1*N, fhat_alpha_bE);
        fhat_alpha_1.intract(2*N, fhat_alpha_aS);
        fhat_alpha_1.intract(3*N, fhat_alpha_bS);

        Vector<double> fhat_alpha_3(2*N);
        fhat_alpha_3.intract(0*N, fhat_alpha_aN);
        fhat_alpha_3.intract(1*N, fhat_alpha_bN);

        Vector<double> fhat_beta_2(4*N);
        fhat_beta_2.intract(0*N, fhat_beta_aW);
        fhat_beta_2.intract(1*N, fhat_beta_bE);
        fhat_beta_2.intract(2*N, fhat_beta_aN);
        fhat_beta_2.intract(3*N, fhat_beta_bN);

        Vector<double> fhat_beta_3(2*N);
        fhat_beta_3.intract(0*N, fhat_beta_aS);
        fhat_beta_3.intract(1*N, fhat_beta_bS);

        // Perform merge operation
        Matrix<double> S = mergeOperationS(T_alpha_33, T_beta_33, T_alpha_31, T_beta_32);
        Matrix<double> T = mergeOperationT(T_alpha_11, T_beta_22, T_alpha_13, T_beta_23, S);
        Vector<double> fhat = mergeOperationF(T_alpha_13, T_beta_23, T_alpha_33, T_beta_33, fhat_alpha_1, fhat_beta_2, fhat_alpha_3, fhat_beta_3);

        // Reorder for square WESN ordering
        Matrix<double> reorder_matrix(8*N, 8*N);
        double z = 0.0;
        Matrix<double> eyeN(N, N, z);
        eyeN.identify();
        reorder_matrix.intract(0*N, 0*N, eyeN);
        reorder_matrix.intract(1*N, 4*N, eyeN);
        reorder_matrix.intract(2*N, 1*N, eyeN);
        reorder_matrix.intract(3*N, 5*N, eyeN);
        reorder_matrix.intract(4*N, 2*N, eyeN);
        reorder_matrix.intract(5*N, 3*N, eyeN);
        reorder_matrix.intract(6*N, 6*N, eyeN);
        reorder_matrix.intract(7*N, 7*N, eyeN);
        T = reorder_matrix * T;
        fhat = reorder_matrix * fhat;

        // @@TODO: Rerorder S matrix (if necessary)

        return merge_values{S, T, fhat};

}

merge_values mergeHorizontal2(
    Matrix<double>& T_alpha,
    Matrix<double>& T_beta,
    Vector<double>& fhat_alpha,
    Vector<double>& fhat_beta,
    int N_levels,
    int level) {

        // Get N from level in tree
        int N = (T_alpha.rows() / 4) * (N_levels - level);

        // Create permutation matrix
        //    Alpha
        Vector<int> pi_alpha({0, 2, 3, 1});
        Matrix<double> P_alpha = blockPermutation<double>(pi_alpha, N);
        Matrix<double> P_alphaT = P_alpha.transpose();

        //    Beta
        Vector<int> pi_beta({1, 2, 3, 0});
        Matrix<double> P_beta = blockPermutation<double>(pi_beta, N);
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

        return merge_values{S, T, fhat};

}

merge_values mergeVertical2(
    Matrix<double>& T_alpha,
    Matrix<double>& T_beta,
    Vector<double>& fhat_alpha,
    Vector<double>& fhat_beta,
    int N_levels,
    int level) {

        // Get N from level in tree
        int N = (T_alpha.rows() / 6) * (N_levels - level - 1);

        // Create permutation matrix
        //    Alpha
        Vector<int> pi_alpha({0, 3, 1, 4, 2, 5});
        Matrix<double> P_alpha = blockPermutation<double>(pi_alpha, N);
        Matrix<double> P_alphaT = P_alpha.transpose();

        //    Beta
        Vector<int> pi_beta({0, 3, 2, 5, 1, 4});
        Matrix<double> P_beta = blockPermutation<double>(pi_beta, N);
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

        // Reorder to return to WESN ordering
        Vector<int> pi_tau({0, 4, 1, 5, 2, 3, 6, 7});
        Matrix<double> P_tau = blockPermutation<double>(pi_tau, N);
        Matrix<double> P_tauT = P_tau.transpose();

        T = P_tau * T;
        T = T * P_tauT;
        fhat = P_tau * fhat;

        return merge_values{S, T, fhat};

}


} // END NAMESPACE hps
