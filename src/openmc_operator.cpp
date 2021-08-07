#include "openmc_operator.h"

OpenMCOperator::OpenMCOperator(int argc, char* argv[], const void* intracomm) {
	// initialize openmc
	openmc_init(argc, argv, intracomm);

	// store the indices of depletable mats for depletion calculation
	// the key of entry is mat id, and the value of entry is index in materials vector
	for (auto& entry : openmc::model::material_map) {
		auto& mat = openmc::model::materials[entry.second];
		if (mat->depletable_) {
			mats_indices_in_depletion_.emplace_back(entry.second);
		}
	}
	std::sort(mats_indices_in_depletion_.begin(), mats_indices_in_depletion_.end());
};

// TODO: consider how to set nuclides in openmc
void OpenMCOperator::InitNuclides(std::vector<std::string>& _nuclides_in_dep_lib) {
	std::set<std::string> nuclides_in_openmc;
	for (const auto& nuclide_name : _nuclides_in_dep_lib) {
		if (openmc::data::nuclide_map.find(nuclide_name) == openmc::data::nuclide_map.end()) continue;
		
		nuclides_in_openmc.insert(nuclide_name);
	}

	// get initial densities using material id 
	for (auto& entry : mats_indices_in_depletion_) {
		// initial nuclides
		int* nuclide_indices;
		double* nuclide_dens;
		const char* nuclide_name;
		int num_nucl;

		openmc_material_get_densities(entry, &nuclide_indices, &nuclide_dens, &num_nucl);


		std::set<std::string> duplicate_nuclides_in_openmc(nuclides_in_openmc);

		for (int i = 0; i < num_nucl; ++i) {
			openmc_nuclide_name(*(nuclide_indices + i), &nuclide_name);
			nuclides_in_openmc.insert(nuclide_name);
			duplicate_nuclides_in_openmc.erase(nuclide_name);
		}

		// give a non-zero value to the nuclides which are not in the initial nuclide list
		std::vector<std::string> nuclide_not_in_openmc(duplicate_nuclides_in_openmc.begin(), duplicate_nuclides_in_openmc.end());
		std::vector<double> nuclide_dens_min(nuclide_not_in_openmc.size());
		for (int i = 0; i < nuclide_dens_min.size(); ++i) {
			nuclide_dens_min[i] = 1.0e-21; 
		}

		auto& mat = openmc::model::materials[entry];
		mat->set_densities(nuclide_not_in_openmc, nuclide_dens_min);
	}

	for (const auto& ele : nuclides_in_openmc) {
		nuclides_in_openmc_.emplace_back(ele);
	}


}

void OpenMCOperator::UnpackTallies() {
	size_t shape[3];
	double* res;
	openmc_tally_results(flux_tally_index_, &res, shape);

	for (int i = 0; i < shape[0]; ++i) {
		res[i];
	}

	for (int j = 0; j < nuclide_in)
};

void OpenMCOperator::InitTallies(std::vector<std::string>& _xs_score) {
	// create a new material filter
	auto material_filter = openmc::Filter::create<openmc::MaterialFilter>();

	// add all cells to the cell filter
	std::vector<int32_t> material_indices;
	for (auto& entry : openmc::model::material_map) {
		material_indices.push_back(entry.second);
	}

	// sort to make sure the cell bins appear in the same
		// order as the test relying on the openmc exe
	std::sort(material_indices.begin(), material_indices.end());
	material_filter->set_materials(material_indices);

	auto energy_filter = openmc::Filter::create<openmc::EnergyFilter>();
	std::vector<double> energy_bins{ 0.0,20.0e6 };
	energy_filter->set_bins(energy_bins);

	// create energy function filter for fission yields tallies
	auto energy_function_filter = openmc::Filter::create<openmc::EnergyFunctionFilter>();
	energy_function_filter->set_data(energy_bins, energy_bins);

	std::vector<openmc::Filter*> filters = { material_filter };

	// create flux tally
	flux_tally_ = openmc::Tally::create();
	flux_tally_->set_filters(filters);
	flux_tally_->set_scores({ "flux" });
	flux_tally_index_ = openmc::model::tallies.size() - 1;

	// create xs tally
	xs_tally_ = openmc::Tally::create();
	xs_tally_->set_filters(filters);
	xs_tally_->set_scores(_xs_score);
	xs_tally_->set_nuclides(nuclides_in_openmc_);
	xs_tally_index_ = openmc::model::tallies.size() - 1;

	// create a heat tally
	heat_tally_ = openmc::Tally::create();
	heat_tally_->set_scores({"heating_local"});
	heat_tally_->set_nuclides(nuclides_in_openmc_);
	heat_tally_index_ = openmc::model::tallies.size() - 1;

	// create a fy tally
	fy_tally_ = openmc::Tally::create();
	filters.emplace_back(energy_filter);
	fy_tally_->set_filters(filters);
	fy_tally_->set_scores({ "fission" });
	fy_tally_->set_nuclides(nuclides_in_openmc_);
	fy_tally_index_ = openmc::model::tallies.size() - 1;

	// create a weighted tally for fy
	weighted_tally_ = openmc::Tally::create();
	filters.emplace_back(energy_function_filter);
	weighted_tally_->set_filters(filters);
	weighted_tally_->set_scores({ "fission" });
	weighted_tally_->set_nuclides(nuclides_in_openmc_);
	weighted_tally_index_ = openmc::model::tallies.size() - 1;

};

void OpenMCOperator::SetNuclideDensities(const int& _mat_id, std::vector<double>& _nuclide_dens) {
	int mat_index = openmc::model::material_map[_mat_id];
	auto& mat = openmc::model::materials[mat_index];

	// the order of _nuclide_dens must be consistent with nuclides_in_openmc_
	mat->set_densities(nuclides_in_openmc_, _nuclide_dens);
}