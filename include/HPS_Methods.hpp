#ifndef HPS_METHODS_HPP_
#define HPS_METHODS_HPP_

#include "HPS_Base.hpp"
#include "HPS_Patch.hpp"
#include "HPS_PatchSolver.hpp"

namespace hps {

struct merge_values {
    Matrix<double> S;
    Matrix<double> T;
    Vector<double> fhat;
};

Matrix<double> mergeOperationS(Matrix<double> T_alpha_33, Matrix<double> T_beta_33, Matrix<double> T_alpha_31, Matrix<double> T_beta_32);
Matrix<double> mergeOperationT(Matrix<double> T_alpha_11, Matrix<double> T_beta_22, Matrix<double> T_alpha_13, Matrix<double> T_beta_23, Matrix<double> S);
// std::pair<Matrix<double>, Matrix<double>> mergeHorizontal(const Matrix<double>& T_level, int N_levels, int level);
// std::pair<Matrix<double>, Matrix<double>> mergeVertical(const Matrix<double>& T_level, int N_levels, int level);
merge_values mergeHorizontal(const Matrix<double>& T_level, int N_levels, int level);
merge_values mergeVertical(const Matrix<double>& T_level, int N_levels, int level);
merge_values mergeHorizontal(const Matrix<double>& T_alpha, const Matrix<double>& T_beta, const Vector<double>& fhat_alpha, const Vector<double>& fhat_beta, int N_levels, int level);
merge_values mergeVertical(const Matrix<double>& T_alpha, const Matrix<double>& T_beta, const Vector<double>& fhat_alpha, const Vector<double>& fhat_beta, int N_levels, int level);

}


#endif // HPS_METHODS_HPP_
