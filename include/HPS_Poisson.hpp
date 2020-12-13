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

class PoissonProblem {

public:

	int ID;
	double Ax;
	double Bx;
	double Ay;
	double By;

	PoissonProblem(int problem_ID, double Ax, double Bx, double Ay, double By);
	double u(double x, double y);
	double f(double x, double y);
	double gEast(double y);
	double gWest(double y);
	double gSouth(double x);
	double gNorth(double x);
	double hEast(double y);
	double hWest(double y);
	double hSouth(double x);
	double hNorth(double x);

};

} // Namespace: hps

#endif // HPS_POISSON_HPP_
