#include "HPS_Poisson.hpp"

namespace hps {

double poisson_u(double x, double y, int PROB_ID) {

	if (PROB_ID == LINEAR) {
		return 1.0;
	}
	else if (PROB_ID == QUAD) {
		return x*x + y*y;
	}
	else if (PROB_ID == POLY) {
		return y*(1 - y)*(x*x*x);
	}
	else if (PROB_ID == SIN) {
		return sin(2*M_PI*(x - 0.5)) * sin(2*M_PI*(y - 0.5));
	}
	else if (PROB_ID == MIX) {
		return (4 + 3*x + 2*y - pow(x,2) + 4*pow(y,2) + 2*x*y + 3*sin(M_PI*x) + 3*cos(M_PI*y) + 2*sin(pow(x,2) + pow(y,2)));
	}
	else {
		printf("ERROR : INVALID PROB_ID\n");
		exit(1);
	}

}

double poisson_f(double x, double y, int PROB_ID) {

	if (PROB_ID == LINEAR) {
		return 0.0;
	}
	else if (PROB_ID == QUAD) {
		return 4.0;
	}
	else if (PROB_ID == POLY) {
		return 6*x*y*(1 - y) - 2*(x*x*x);
	}
	else if (PROB_ID == SIN) {
		return (-8.0*M_PI*M_PI)*sin(2*M_PI*(x - 0.5)) * sin(2*M_PI*(y - 0.5));
	}
	else if (PROB_ID == MIX) {
		return (6 - 3*pow(M_PI,2)*cos(M_PI*y) + 8*cos(pow(x,2) + pow(y,2)) - 3*pow(M_PI,2)*sin(M_PI*x) - 8*pow(x,2)*sin(pow(x,2) + pow(y,2)) - 8*pow(y,2)*sin(pow(x,2) + pow(y,2)));
	}
	else {
		printf("ERROR : INVALID PROB_ID\n");
		exit(1);
	}

}

double poisson_g(double x, double y, int PROB_ID, int BC_SIDE) {

	if (PROB_ID == LINEAR) {
		return 1.0;
	}
	else if (PROB_ID == QUAD) {
		if (BC_SIDE == WEST) {
			return y*y;
		}
		else if (BC_SIDE == EAST) {
			return 1 + y*y;
		}
		else if (BC_SIDE == SOUTH) {
			return x*x;
		}
		else if (BC_SIDE == NORTH) {
			return 1 + x*x;
		}
		else {
			printf("ERROR : INVALID BC_SIDE\n");
			exit(1);
		}
	}
	else if (PROB_ID == POLY) {
		if (BC_SIDE == EAST) {
			return y*(1 - y);
		}
		else if (BC_SIDE == WEST || BC_SIDE == SOUTH || BC_SIDE == NORTH) {
			return 0.0;
		}
		else {
			printf("ERROR : INVALID BC_SIDE\n");
			exit(1);
		}
	}
	else if (PROB_ID == SIN) {
		return 0.0;
	}
	else if (PROB_ID == MIX) {
		if (BC_SIDE == WEST) {
			return (4 + 2*y + 4*pow(y,2) + 3*cos(M_PI*y) + 2*sin(pow(y,2)));
		}
		else if (BC_SIDE == EAST) {
			return (6 + 4*y + 4*pow(y,2) + 3*cos(M_PI*y) + 2*sin(1 + pow(y,2)));
		}
		else if (BC_SIDE == SOUTH) {
			return (7 + 3*x - pow(x,2) + 3*sin(M_PI*x) + 2*sin(pow(x,2)));
		}
		else if (BC_SIDE == NORTH) {
			return (7 + 5*x - pow(x,2) + 3*sin(M_PI*x) + 2*sin(1 + pow(x,2)));
		}
		else {
			printf("ERROR : INVALID BC_SIDE\n");
			exit(1);
		}
	}
	else {
		printf("ERROR : INVALID PROB_ID\n");
		exit(1);
	}

}

double poisson_h(double x, double y, int PROB_ID, int BC_SIDE) { // Change to return gradient; dot with normal vector

	if (PROB_ID == LINEAR) {
		return 0.0;
	}
	else if (PROB_ID == QUAD) {
		if (BC_SIDE  == WEST || BC_SIDE == SOUTH) {
			return 0.0;
		}
		else {
			return 2.0;
		}
	}
	else if (PROB_ID == POLY) {
		if (BC_SIDE == WEST) {
			return 0.0;
		}
		else if (BC_SIDE == EAST) {
			return 3*y*(1 - y);
		}
		else if (BC_SIDE == SOUTH) {
			return y*y*y;
		}
		else if (BC_SIDE == NORTH) {
			return -(y*y*y);
		}
	}
	else if (PROB_ID == SIN) {
		if (BC_SIDE == WEST || BC_SIDE == EAST) {
			return -2*M_PI*sin(2*M_PI*(y - 0.5));
		}
		else if (BC_SIDE == SOUTH || BC_SIDE == NORTH) {
			return -2*M_PI*sin(2*M_PI*(x - 0.5));
		}
	}
	else if (PROB_ID == MIX) {
		if (BC_SIDE == WEST) {
			return (3 + 3*M_PI + 2*y);
		}
		else if (BC_SIDE == EAST) {
			return (1 - 3*M_PI + 2*y + 4*cos(1 + pow(y,2)));
		}
		else if (BC_SIDE == SOUTH) {
			return (2 + 2*x);
		}
		else if (BC_SIDE == NORTH) {
			return (10 + 2*x + 4*cos(1 + pow(x,2)));
		}
		else {
			printf("ERROR : INVALID BC_SIDE\n");
			exit(1);
		}
	}
	else {
		printf("ERROR : INVALID PROB_ID\n");
		exit(1);
	}

}

PoissonProblem::PoissonProblem(int problem_ID, double Ax, double Bx, double Ay, double By) : ID(problem_ID), Ax(Ax), Bx(Bx), Ay(Ay), By(By) {}

double PoissonProblem::u(double x, double y) {
	return 2*x*y + pow(x,2) + pow(y,2);
}

double PoissonProblem::f(double x, double y) {
	return 4.0;
}

double PoissonProblem::gEast(double y) {
	return pow(Ax,2) + 2*Ax*y + pow(y,2);
}

double PoissonProblem::gWest(double y) {
	return pow(Bx,2) + 2*Bx*y + pow(y,2);
}

double PoissonProblem::gSouth(double x) {
	return pow(Ay,2) + 2*Ay*x + pow(x,2);
}

double PoissonProblem::gNorth(double x) {
	return pow(By,2) + 2*By*x + pow(x,2);
}

double PoissonProblem::hEast(double y) {
	return 2*Ax + 2*y;
}

double PoissonProblem::hWest(double y) {
	return 2*Bx + 2*y;
}

double PoissonProblem::hSouth(double x) {
	return 2*Ay + 2*x;
}

double PoissonProblem::hNorth(double x) {
	return 2*By + 2*x;
}

}
