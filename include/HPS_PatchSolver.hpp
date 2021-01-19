#ifndef HPS_PATCHSOLVER_HPP_
#define HPS_PATCHSOLVER_HPP_

#include "HPS_Base.hpp"
#include <iostream>

extern "C" {
	void hstcrt_(double* A, double* B, int* M, int* MBDCND, double* BDA, double* BDB, double* C, double* D, int* N, int* NBDCND, double* BDC, double* BDD, double* ELMBDA, double* F, int* IDIMF, double* PERTRB, int* IERROR, double* W);
}

namespace hps {

hps::Matrix<double> mapSolution(hps::CellGrid<double, 2>& grid, hps::Vector<double>& g, hps::Matrix<double>& f, double lambda);
hps::Vector<double> mapDirichletToNeumann(hps::CellGrid<double, 2> grid, hps::Vector<double> g, hps::Matrix<double> f, double lambda);
hps::Matrix<double> buildDirichletToNeumann(hps::CellGrid<double, 2> grid, double lambda);

} // Namespace hps

#endif // HPS_PATCHSOLVER_HPP_
