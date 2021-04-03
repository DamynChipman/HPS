#ifndef HPS_MERGE_HPP_
#define HPS_MERGE_HPP_

#include "HPS_Base.hpp"
#include "HPS_Patch.hpp"
#include "HPS_PatchSolver.hpp"

namespace hps {

/**
 * Struct containing patch data
 * @data S     Solution operator matrix
 * @data T     DtN operator matrix
 * @data fhat  fhat vector
 * @data w     w vector
 */
struct merge_values {
    Matrix<double> S;
    Matrix<double> T;
    Vector<double> fhat;
    Vector<double> w;
};

/**
 * Performs the merge operation (linear algebra) for the solution operator S.
 * @NOTE: This HAS been unit tested against Python and passed
 * @param  T_alpha_33 Matrix for T_alpha_33
 * @param  T_beta_33  Matrix for T_beta_33
 * @param  T_alpha_31 Matrix for T_alpha_31
 * @param  T_beta_32  Matrix for T_beta_32
 * @return            Computed solution operator S
 */
Matrix<double> mergeOperationS(Matrix<double> T_alpha_33, Matrix<double> T_beta_33, Matrix<double> T_alpha_31, Matrix<double> T_beta_32);

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
Matrix<double> mergeOperationT(Matrix<double> T_alpha_11, Matrix<double> T_beta_22, Matrix<double> T_alpha_13, Matrix<double> T_beta_23, Matrix<double> S);

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
Vector<double> mergeOperationF(Matrix<double> T_alpha_13, Matrix<double> T_beta_23, Matrix<double> T_alpha_33, Matrix<double> T_beta_33, Vector<double> fhat_alpha_1, Vector<double> fhat_beta_2, Vector<double> fhat_alpha_3, Vector<double> fhat_beta_3);

/**
 * Performs the merge operation (linear algebra) for the w vector
 * @param  T_alpha_33   Matrix for T_alpha_33
 * @param  T_beta_33    Matrix for T_beta_33
 * @param  fhat_alpha_3 Vector for fhat_alpha_3
 * @param  fhat_beta_3  Vector for fhat_beta_3
 * @return              Computed w vector
 */
Vector<double> mergeOperationW(Matrix<double> T_alpha_33, Matrix<double> T_beta_33, Vector<double> fhat_alpha_3, Vector<double> fhat_beta_3);

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
merge_values mergeHorizontal(Matrix<double>& T_alpha, Matrix<double>& T_beta, Vector<double>& fhat_alpha, Vector<double>& fhat_beta, int N_levels, int level);

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
merge_values mergeVertical(Matrix<double>& T_alpha, Matrix<double>& T_beta, Vector<double>& fhat_alpha, Vector<double>& fhat_beta, int N_levels, int level);

}


#endif // HPS_MERGE_HPP_
