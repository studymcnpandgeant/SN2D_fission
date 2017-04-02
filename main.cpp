#include "Solver.h"


#include <iostream>

int main() {

	Quadrature* quad = new Quadrature(4);
//	quad->initializeArray();
	Solver* solver = new Solver(100,100,4);
	solver->setQuadrature(quad);
	solver->startIteration();

	std::cout << "Residual = " << solver->getResidual() << std::endl;

	return 0;

}
