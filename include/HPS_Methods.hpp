#ifndef HPS_METHODS_HPP_
#define HPS_METHODS_HPP_

#include "HPS_Base.hpp"
#include "HPS_Patch.hpp"
#include "HPS_PatchSolver.hpp"

namespace hps {

hps::Matrix<double> mergeOperationS(hps::Matrix<double> T_alpha_33, hps::Matrix<double> T_beta_33, hps::Matrix<double> T_alpha_31, hps::Matrix<double> T_beta_32);
hps::Matrix<double> mergeOperationT(hps::Matrix<double> T_alpha_11, hps::Matrix<double> T_beta_22, hps::Matrix<double> T_alpha_13, hps::Matrix<double> T_beta_23, hps::Matrix<double> S);
std::pair<hps::Matrix<double>, hps::Matrix<double>> mergeHorizontal(const hps::Matrix<double>& T_level, int N_levels, int level);
std::pair<hps::Matrix<double>, hps::Matrix<double>> mergeVertical(const hps::Matrix<double>& T_level, int N_levels, int level);

}


#endif // HPS_METHODS_HPP_
