//#include "nuclide_system.h"
#include "modec_solver.h"

ModecSolver& ModecSolver::getModecSolver()
{
    static ModecSolver instance;
    return instance;
};

void ModecSolver::SetSolverParams(const std::string& _name, const int& _major_order, const int& _minor_order) {
    std::string solver_name_ = _name;
    // convert solver name to upper case
    std::transform(solver_name_.begin(), solver_name_.end(), solver_name_.begin(), ::toupper);
    
    auto iter = MODEC::SolverIndicator.find(solver_name_);
    if (iter == MODEC::SolverIndicator.end()) {
        spdlog::warn("There is no existed solver named {} !!\n IpfCRAM48 is used.", _name);
    }
    else {
        solver_id_ = iter->second;
    }

    major_order_ = _major_order;
    minor_order_ = _minor_order;
};

void ModecSolver::Execute(SparseMatrix<double>& _burnup_matrix, SymbolicSpMat& _sym_spmat, double& _time, std::vector<double>& _nuclide_densities) {
    switch (solver_id_) {
    // using PfdCRAM solver
    case 1:
        switch (major_order_) {
        case 14:
            PfdCRAM<14>::Solve(_burnup_matrix, _sym_spmat, _time, _nuclide_densities);
            break;
        case 16:
            PfdCRAM<16>::Solve(_burnup_matrix, _sym_spmat, _time, _nuclide_densities);
            break;
        default:
            spdlog::warn("No existed {}th-order PFD-CRAM method!! 16th-order PFD-CRAM is used instead !! \n", major_order_);
            PfdCRAM<16>::Solve(_burnup_matrix, _sym_spmat, _time, _nuclide_densities);
            break;
        }
        break;

    // using IpfCRAM solver
    case 2:
        switch (major_order_) {
        case 32:
            IpfCRAM<32>::Solve(_burnup_matrix, _sym_spmat, _time, _nuclide_densities);
            break;
        case 48:
            IpfCRAM<48>::Solve(_burnup_matrix, _sym_spmat, _time, _nuclide_densities);
            break;
        default:
            spdlog::warn("No existed {}th-order Ipf-CRAM method!! 48th-order Ipf-CRAM is used instead !! \n", major_order_);
            IpfCRAM<48>::Solve(_burnup_matrix, _sym_spmat, _time, _nuclide_densities);
            break;
        }
        break;

    // using PfdPRAM solver
    case 3:
        switch (major_order_) {
        case 16:
            PfdPRAM<16,4>::Solve(_burnup_matrix, _sym_spmat, _time, _nuclide_densities);
            break;
        case 32:
            PfdPRAM<32,8>::Solve(_burnup_matrix, _sym_spmat, _time, _nuclide_densities);
            break;
        default:
            spdlog::warn("No existed {}th-order PFD-PRAM method!! 32th-order PFD-PRAM is used instead !! \n", major_order_);
            PfdPRAM<32, 8>::Solve(_burnup_matrix, _sym_spmat, _time, _nuclide_densities);
            break;
        }
        break;

    // using IpfPRAM solver
    case 4:
        switch (major_order_) {
        case 16:
            IpfPRAM<16, 4>::Solve(_burnup_matrix, _sym_spmat, _time, _nuclide_densities);
            break;
        case 32:
            IpfPRAM<32, 8>::Solve(_burnup_matrix, _sym_spmat, _time, _nuclide_densities);
            break;
        case 48:
            IpfPRAM<48, 16>::Solve(_burnup_matrix, _sym_spmat, _time, _nuclide_densities);
            break;
        case 64:
            IpfPRAM<64, 32>::Solve(_burnup_matrix, _sym_spmat, _time, _nuclide_densities);
            break;
        default:
            spdlog::warn("No existed {}th-order Ipf-PRAM method!! 64th-order Ipf-PRAM is used instead !! \n", major_order_);
            IpfPRAM<64, 32>::Solve(_burnup_matrix, _sym_spmat, _time, _nuclide_densities);
            break;
        }
        break;
    }
};
