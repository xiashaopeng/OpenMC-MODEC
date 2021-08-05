#include "modec_solver.h"

int main() {
	SparseMatrix<double> _burnup_matrix;
	SymbolicSpMat _sym_spmat;
	double _time;
	std::vector<double> _nuclide_densities;

	auto& a = ModecSolver::getModecSolver();
	a.SetSolverParams("PfdCRAM", 16);
	a.Execute(_burnup_matrix, _sym_spmat, _time, _nuclide_densities);
	return 0;
}