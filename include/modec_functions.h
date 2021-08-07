#ifndef MODEC_FUNCTIONS_H_
#define MODEC_FUNCTIONS_H_

#include "dep_zones.h"
#include "modec_solver.h"
#include "openmc_operator.h"

#include "tinyxml2.h"
#include "spdlog.h"

namespace MODEC {

// initialize modec
void InitModec(const std::string& _lib_filename);

// create depletion zones
// TODO: set initial densities
void InitDepZones(DepZones& _dep_zones, std::vector<int>& _zone_indices);

// get burnup matrix A
SparseMatrix<double> AssembleMat(NuclideSystem& _nuclide_system, const double& norm_factor);

// calculation phi-function 
// \varphi_0(tA) \ctimes \vec{b}_0 + \varphi_1(tA) \ctimes t\vec{b}_1 + \varphi_2(tA) \ctimes t^2\vec{b}_2 + \varphi_3(tA) \ctimes t^3\vec{b}_3 + ...
void CalcPhiFunction(SparseMatrix<double> _burnup_matrix, SymbolicSpMat& sym_spmat, std::map<int, std::vector<double> >& _vec_b);


};

#endif