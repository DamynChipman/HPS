#include "HPS_Poisson.hpp"

namespace hps {

PoissonProblem::PoissonProblem(int problem_ID, double Ax, double Bx, double Ay, double By) : ID(problem_ID), Ax(Ax), Bx(Bx), Ay(Ay), By(By) {}

double PoissonProblem::u(double x, double y) {
	switch(ID) {
		case CONSTANT:
			return 1.0;
		case LINEAR:
			return x + y;
		case LAPLACE:
			return y*sin(2*M_PI*x) + x*cos(2*M_PI*y) + 4;
		case QUAD:
			return x*x + y*y + 2*x*y;
		case POLY:
			return y*y*x*x*x*x;
		case TRIG:
			return sin(2*M_PI*x) * sin(2*M_PI*y);
		default:
			throw std::invalid_argument("[PoissonProblem::u] Invalid problem_ID.");
	}
}

double PoissonProblem::f(double x, double y) {
	switch(ID) {
		case CONSTANT:
			return 0.0;
		case LINEAR:
			return 0.0;
		case LAPLACE:
			return 0.0;
		case QUAD:
			return 4.0;
		case POLY:
			return 2.0*x*x*(6.0*y*y + x*x);
		case TRIG:
			return -2.0*pow(2.0*M_PI,2)*u(x,y);
		default:
			throw std::invalid_argument("[PoissonProblem::u] Invalid problem_ID.");
	}

}

double PoissonProblem::dudx(double x, double y) {
	switch(ID) {
		case CONSTANT:
			return 0.0;
		case LINEAR:
			return 1.0;
		case QUAD:
			return 2.0*x + 2.0*y;
		case POLY:
			return 4.0*y*y*x*x*x;
		case TRIG:
			return 2.0*M_PI*cos(2.0*M_PI*x) * sin(2.0*M_PI*y);
		default:
			throw std::invalid_argument("[PoissonProblem::u] Invalid problem_ID.");
	}
}

double PoissonProblem::dudy(double x, double y) {
	switch(ID) {
		case CONSTANT:
			return 0.0;
		case LINEAR:
			return 1.0;
		case QUAD:
			return 2.0*y + 2*x;
		case POLY:
			return 2.0*y*x*x*x*x;
		case TRIG:
			return sin(2*M_PI*x) * 2.0*M_PI*cos(2*M_PI*y);
		default:
			throw std::invalid_argument("[PoissonProblem::u] Invalid problem_ID.");
	}
}

// double PoissonProblem::gEast(double y) {
// 	return pow(Ax,2) + 2*Ax*y + pow(y,2);
// }
//
// double PoissonProblem::gWest(double y) {
// 	return pow(Bx,2) + 2*Bx*y + pow(y,2);
// }
//
// double PoissonProblem::gSouth(double x) {
// 	return pow(Ay,2) + 2*Ay*x + pow(x,2);
// }
//
// double PoissonProblem::gNorth(double x) {
// 	return pow(By,2) + 2*By*x + pow(x,2);
// }
//
// double PoissonProblem::hEast(double y) {
// 	return 2*Ax + 2*y;
// }
//
// double PoissonProblem::hWest(double y) {
// 	return 2*Bx + 2*y;
// }
//
// double PoissonProblem::hSouth(double x) {
// 	return 2*Ay + 2*x;
// }
//
// double PoissonProblem::hNorth(double x) {
// 	return 2*By + 2*x;
// }

}
