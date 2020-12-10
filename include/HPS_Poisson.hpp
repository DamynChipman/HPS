#ifndef HPS_POISSON_HPP_
#define HPS_POISSON_HPP_

#include <iostream>
#include <cmath>

namespace hps {

// ===== Enums =====
enum PROBLEM_TYPE {
	LINEAR,
	QUAD,
	POLY,
	SIN,
	MIX
};

enum BOUNDARY_SIDE {
	WEST,
	EAST,
	SOUTH,
	NORTH
};

// ===== Function Definitions =====
double poisson_u(double x, double y, int PROB_ID);
double poisson_f(double x, double y, int PROB_ID);
double poisson_g(double x, double y, int PROB_ID, int BC_SIDE);
double poisson_h(double x, double y, int PROB_ID, int BC_SIDE);
    
} // Namespace: hps

#endif // HPS_POISSON_HPP_
