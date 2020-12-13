#include "HPS_Methods.hpp"

namespace hps {

/**
 * Performs the merge operation (linear algebra) for the solution operator S.
 * @param  T_alpha_33 Matrix for T_alpha_33
 * @param  T_beta_33  Matrix for T_beta_33
 * @param  T_alpha_31 Matrix for T_alpha_31
 * @param  T_beta_32  Matrix for T_beta_32
 * @return            Computed solution operator S
 */
hps::Matrix<double> mergeOperationS(hps::Matrix<double> T_alpha_33, hps::Matrix<double> T_beta_33, hps::Matrix<double> T_alpha_31, hps::Matrix<double> T_beta_32) {

    hps::Matrix<double> S_RHS(T_alpha_31.rows(), T_alpha_31.cols() + T_beta_32.cols());
    T_alpha_31.negate();
    S_RHS.intract(0, 0, T_alpha_31);
    S_RHS.intract(0, T_alpha_31.cols(), T_beta_32);
    hps::Matrix<double> S = hps::solve(T_alpha_33 - T_beta_33, S_RHS);
    return S;

}

/**
 * Performs the merge operation (linear algebra) for the DtN operator T.
 * @param  T_alpha_11 Matrix for T_alpha_11
 * @param  T_beta_22  Matrix for T_beta_22
 * @param  T_alpha_13 Matrix for T_alpha_13
 * @param  T_beta_23  Matrix for T_beta_23
 * @param  S          Matrix for solution operator
 * @return            Computed DtN operator T
 */
hps::Matrix<double> mergeOperationT(hps::Matrix<double> T_alpha_11, hps::Matrix<double> T_beta_22, hps::Matrix<double> T_alpha_13, hps::Matrix<double> T_beta_23, hps::Matrix<double> S) {

    hps::Matrix<double> T(T_alpha_11.rows() + T_beta_22.rows(), T_alpha_11.cols() + T_beta_22.cols());
    T.intract(0, 0, T_alpha_11);
    T.intract(T_alpha_11.rows(), T_alpha_11.cols(), T_beta_22);
    hps::Matrix<double> T_RHS(T_alpha_13.rows() + T_beta_23.rows(), T_alpha_13.cols());
    T_RHS.intract(0, 0, T_alpha_13);
    T_RHS.intract(T_alpha_13.rows(), 0, T_beta_23);
    T += T_RHS*S;
    return T;

}

/**
 * Performs a horizontal merge for a provided DtN operator `T_level` at the specified `level`. Extracts the blocks corresponding to each side of the patch, intracts them into the merge operation matrices, and calls `mergeOperationS` and `mergeOperationT` to do the linear algebra.
 * @param T_level  Matrix for DtN operator at the specified level. Should be a (4*N x 4*N) matrix for each level away from the leaf level
 * @param N_levels Number of levels in the tree
 * @param level    Current level. Should be even for the horizontal merge
 * @return         Returns a std::pair of matrices (S,T)
 */
std::pair<hps::Matrix<double>, hps::Matrix<double>> mergeHorizontal(const hps::Matrix<double>& T_level, int N_levels, int level) {

    // Get number of entries in each section
    int N = (T_level.rows() / 4) * (N_levels - level);

    // Unpack input matrix into blocks
    hps::Matrix<double> T_WW = T_level.extract(0*N, 0*N, N, N);
    hps::Matrix<double> T_WE = T_level.extract(0*N, 1*N, N, N);
    hps::Matrix<double> T_WS = T_level.extract(0*N, 2*N, N, N);
    hps::Matrix<double> T_WN = T_level.extract(0*N, 3*N, N, N);

    hps::Matrix<double> T_EW = T_level.extract(1*N, 0*N, N, N);
    hps::Matrix<double> T_EE = T_level.extract(1*N, 1*N, N, N);
    hps::Matrix<double> T_ES = T_level.extract(1*N, 2*N, N, N);
    hps::Matrix<double> T_EN = T_level.extract(1*N, 3*N, N, N);

    hps::Matrix<double> T_SW = T_level.extract(2*N, 0*N, N, N);
    hps::Matrix<double> T_SE = T_level.extract(2*N, 1*N, N, N);
    hps::Matrix<double> T_SS = T_level.extract(2*N, 2*N, N, N);
    hps::Matrix<double> T_SN = T_level.extract(2*N, 3*N, N, N);

    hps::Matrix<double> T_NW = T_level.extract(3*N, 0*N, N, N);
    hps::Matrix<double> T_NE = T_level.extract(3*N, 1*N, N, N);
    hps::Matrix<double> T_NS = T_level.extract(3*N, 2*N, N, N);
    hps::Matrix<double> T_NN = T_level.extract(3*N, 3*N, N, N);

    // Create partitioned matrices
    hps::Matrix<double> T_alpha_11(3*N, 3*N);
    T_alpha_11.intract(0*N, 0*N, T_WW);
    T_alpha_11.intract(0*N, 1*N, T_WS);
    T_alpha_11.intract(0*N, 2*N, T_WN);
    T_alpha_11.intract(1*N, 0*N, T_SW);
    T_alpha_11.intract(1*N, 1*N, T_SS);
    T_alpha_11.intract(1*N, 2*N, T_SN);
    T_alpha_11.intract(2*N, 0*N, T_NW);
    T_alpha_11.intract(2*N, 1*N, T_NS);
    T_alpha_11.intract(2*N, 2*N, T_NN);

    hps::Matrix<double> T_alpha_13(3*N, 1*N);
    T_alpha_13.intract(0*N, 0*N, T_WE);
    T_alpha_13.intract(1*N, 0*N, T_SE);
    T_alpha_13.intract(2*N, 0*N, T_NE);

    hps::Matrix<double> T_alpha_31(1*N, 3*N);
    T_alpha_31.intract(0*N, 0*N, T_EW);
    T_alpha_31.intract(0*N, 1*N, T_ES);
    T_alpha_31.intract(0*N, 2*N, T_EN);

    hps::Matrix<double> T_alpha_33(1*N, 1*N);
    T_alpha_33.intract(0*N, 0*N, T_EE);

    hps::Matrix<double> T_beta_22(3*N, 3*N);
    T_beta_22.intract(0*N, 0*N, T_EE);
    T_beta_22.intract(0*N, 1*N, T_ES);
    T_beta_22.intract(0*N, 2*N, T_EN);
    T_beta_22.intract(1*N, 0*N, T_SE);
    T_beta_22.intract(1*N, 1*N, T_SS);
    T_beta_22.intract(1*N, 2*N, T_SN);
    T_beta_22.intract(2*N, 0*N, T_NE);
    T_beta_22.intract(2*N, 1*N, T_NS);
    T_beta_22.intract(2*N, 2*N, T_NN);

    hps::Matrix<double> T_beta_23(3*N, 1*N);
    T_beta_23.intract(0*N, 0*N, T_EW);
    T_beta_23.intract(1*N, 0*N, T_SW);
    T_beta_23.intract(2*N, 0*N, T_NW);

    hps::Matrix<double> T_beta_32(1*N, 3*N);
    T_beta_32.intract(0*N, 0*N, T_WE);
    T_beta_32.intract(0*N, 1*N, T_WS);
    T_beta_32.intract(0*N, 2*N, T_WN);

    hps::Matrix<double> T_beta_33(1*N, 1*N);
    T_beta_33.intract(0*N, 0*N, T_WW);

    // Perform merge operation
    hps::Matrix<double> S = mergeOperationS(T_alpha_33, T_beta_33, T_alpha_31, T_beta_32);
    hps::Matrix<double> T = mergeOperationT(T_alpha_11, T_beta_22, T_alpha_13, T_beta_23, S);

    return std::pair<hps::Matrix<double>, hps::Matrix<double>>(S, T);
}

/**
 * Performs a vertical merge for a provided DtN operator `T_level` at the specified `level`. Extracts the blocks corresponding to each side of the patch, intracts them into the merge operation matrices, and calls `mergeOperationS` and `mergeOperationT` to do the linear algebra. After the linear algebra, reorders S and T back into the WESN ordering for a square (via a permutation matrix multiplication).
 * @param T_level  Matrix for DtN operator at the specified level. Should be a (6*N x 6*N) matrix for each level away from the leaf level
 * @param N_levels Number of levels in the tree
 * @param level    Current level. Should be odd for the vertical merge
 * @return         Returns a std::pair of matrices (S,T)
 */
std::pair<hps::Matrix<double>, hps::Matrix<double>> mergeVertical(const hps::Matrix<double>& T_level, int N_levels, int level) {

    // Get number of entries in each section
    int N = (T_level.rows() / 6) * (N_levels - level - 1);

    // Unpack input matrix into blocks
    hps::Matrix<double> T_WW_aa = T_level.extract(0*N, 0*N, N, N);
    hps::Matrix<double> T_WS_aa = T_level.extract(0*N, 1*N, N, N);
    hps::Matrix<double> T_WN_aa = T_level.extract(0*N, 2*N, N, N);
    hps::Matrix<double> T_WE_ab = T_level.extract(0*N, 3*N, N, N);
    hps::Matrix<double> T_WS_ab = T_level.extract(0*N, 4*N, N, N);
    hps::Matrix<double> T_WN_ab = T_level.extract(0*N, 5*N, N, N);

    hps::Matrix<double> T_SW_aa = T_level.extract(1*N, 0*N, N, N);
    hps::Matrix<double> T_SS_aa = T_level.extract(1*N, 1*N, N, N);
    hps::Matrix<double> T_SN_aa = T_level.extract(1*N, 2*N, N, N);
    hps::Matrix<double> T_SE_ab = T_level.extract(1*N, 3*N, N, N);
    hps::Matrix<double> T_SS_ab = T_level.extract(1*N, 4*N, N, N);
    hps::Matrix<double> T_SN_ab = T_level.extract(1*N, 5*N, N, N);

    hps::Matrix<double> T_NW_aa = T_level.extract(2*N, 0*N, N, N);
    hps::Matrix<double> T_NS_aa = T_level.extract(2*N, 1*N, N, N);
    hps::Matrix<double> T_NN_aa = T_level.extract(2*N, 2*N, N, N);
    hps::Matrix<double> T_NE_ab = T_level.extract(2*N, 3*N, N, N);
    hps::Matrix<double> T_NS_ab = T_level.extract(2*N, 4*N, N, N);
    hps::Matrix<double> T_NN_ab = T_level.extract(2*N, 5*N, N, N);

    hps::Matrix<double> T_EW_ba = T_level.extract(3*N, 0*N, N, N);
    hps::Matrix<double> T_ES_ba = T_level.extract(3*N, 1*N, N, N);
    hps::Matrix<double> T_EN_ba = T_level.extract(3*N, 2*N, N, N);
    hps::Matrix<double> T_EE_bb = T_level.extract(3*N, 3*N, N, N);
    hps::Matrix<double> T_ES_bb = T_level.extract(3*N, 4*N, N, N);
    hps::Matrix<double> T_EN_bb = T_level.extract(3*N, 5*N, N, N);

    hps::Matrix<double> T_SW_ba = T_level.extract(4*N, 0*N, N, N);
    hps::Matrix<double> T_SS_ba = T_level.extract(4*N, 1*N, N, N);
    hps::Matrix<double> T_SN_ba = T_level.extract(4*N, 2*N, N, N);
    hps::Matrix<double> T_SE_bb = T_level.extract(4*N, 3*N, N, N);
    hps::Matrix<double> T_SS_bb = T_level.extract(4*N, 4*N, N, N);
    hps::Matrix<double> T_SN_bb = T_level.extract(4*N, 5*N, N, N);

    hps::Matrix<double> T_NW_ba = T_level.extract(5*N, 0*N, N, N);
    hps::Matrix<double> T_NS_ba = T_level.extract(5*N, 1*N, N, N);
    hps::Matrix<double> T_NN_ba = T_level.extract(5*N, 2*N, N, N);
    hps::Matrix<double> T_NE_bb = T_level.extract(5*N, 3*N, N, N);
    hps::Matrix<double> T_NS_bb = T_level.extract(5*N, 4*N, N, N);
    hps::Matrix<double> T_NN_bb = T_level.extract(5*N, 5*N, N, N);

    // Create partitioned matricies
    hps::Matrix<double> T_alpha_11(4*N, 4*N);
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

    hps::Matrix<double> T_alpha_13(4*N, 2*N);
    T_alpha_13.intract(0*N, 0*N, T_WN_aa);
    T_alpha_13.intract(0*N, 1*N, T_WN_ab);
    T_alpha_13.intract(1*N, 0*N, T_EN_ba);
    T_alpha_13.intract(1*N, 1*N, T_EN_bb);
    T_alpha_13.intract(2*N, 0*N, T_SN_aa);
    T_alpha_13.intract(2*N, 1*N, T_SN_ab);
    T_alpha_13.intract(3*N, 0*N, T_SN_ba);
    T_alpha_13.intract(3*N, 1*N, T_SN_bb);

    hps::Matrix<double> T_alpha_31(2*N, 4*N);
    T_alpha_31.intract(0*N, 0*N, T_NW_aa);
    T_alpha_31.intract(0*N, 1*N, T_NE_ab);
    T_alpha_31.intract(0*N, 2*N, T_NS_aa);
    T_alpha_31.intract(0*N, 3*N, T_NS_ab);
    T_alpha_31.intract(1*N, 0*N, T_NW_ba);
    T_alpha_31.intract(1*N, 1*N, T_NE_bb);
    T_alpha_31.intract(1*N, 2*N, T_NS_ba);
    T_alpha_31.intract(1*N, 3*N, T_NS_bb);

    hps::Matrix<double> T_alpha_33(2*N, 2*N);
    T_alpha_33.intract(0*N, 0*N, T_NN_aa);
    T_alpha_33.intract(0*N, 1*N, T_NN_ab);
    T_alpha_33.intract(1*N, 0*N, T_NN_ba);
    T_alpha_33.intract(1*N, 1*N, T_NN_bb);

    hps::Matrix<double> T_beta_22(4*N, 4*N);
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

    hps::Matrix<double> T_beta_23(4*N, 2*N);
    T_beta_23.intract(0*N, 0*N, T_WS_aa);
    T_beta_23.intract(0*N, 1*N, T_WS_ab);
    T_beta_23.intract(1*N, 0*N, T_ES_ba);
    T_beta_23.intract(1*N, 1*N, T_ES_bb);
    T_beta_23.intract(2*N, 0*N, T_NS_aa);
    T_beta_23.intract(2*N, 1*N, T_NS_ab);
    T_beta_23.intract(3*N, 0*N, T_NS_ba);
    T_beta_23.intract(3*N, 1*N, T_NS_bb);

    hps::Matrix<double> T_beta_32(2*N, 4*N);
    T_beta_32.intract(0*N, 0*N, T_SW_aa);
    T_beta_32.intract(0*N, 1*N, T_SE_ab);
    T_beta_32.intract(0*N, 2*N, T_SN_aa);
    T_beta_32.intract(0*N, 3*N, T_SN_ab);
    T_beta_32.intract(1*N, 0*N, T_SW_ba);
    T_beta_32.intract(1*N, 1*N, T_SE_bb);
    T_beta_32.intract(1*N, 2*N, T_SN_ba);
    T_beta_32.intract(1*N, 3*N, T_SN_bb);

    hps::Matrix<double> T_beta_33(2*N, 2*N);
    T_beta_33.intract(0*N, 0*N, T_SS_aa);
    T_beta_33.intract(0*N, 1*N, T_SS_ab);
    T_beta_33.intract(1*N, 0*N, T_SS_ba);
    T_beta_33.intract(1*N, 1*N, T_SS_bb);

    // Perform merge operation
    hps::Matrix<double> S = mergeOperationS(T_alpha_33, T_beta_33, T_alpha_31, T_beta_32);
    hps::Matrix<double> T = mergeOperationT(T_alpha_11, T_beta_22, T_alpha_13, T_beta_23, S);

    // Reorder for square WESN ordering
    hps::Matrix<double> reorder_matrix(8*N, 8*N);
    double z = 0.0;
    hps::Matrix<double> eyeN(N, N, z);
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


    return std::pair<hps::Matrix<double>, hps::Matrix<double>>(S, T);
}


}
