#ifndef OPENMC_OPERATOR_H_
#define OPENMC_OPERATOR_H_

#include <openmc/simulation.h>

#include "openmc/capi.h"
#include "openmc/nuclide.h"
#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/material.h"
#include "openmc/message_passing.h"
#include "openmc/summary.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_cell.h"
#include "openmc/tallies/filter_energy.h"
#include "openmc/tallies/filter_energyfunc.h"
#include "openmc/tallies/filter_material.h"
#include "openmc/tallies/filter_mu.h"
#include "openmc/tallies/tally.h"

#include "transport_operator.h"
#include "dep_zones.h"

class OpenMCOperator : public TransportOperator {
public:
	OpenMCOperator(int argc, char* argv[], const void* intracomm);

	void InitializeBatchSimulation() {
		// initialize openmc simulation
		openmc_simulation_init();
	}

	void SimulateBatches(const int& num_batches) {
		// simulate batches
		int status;
		for (int i = 0; i < num_batches; ++i) {
			openmc_next_batch(&status);
		}
	}

	void InitNuclides(std::vector<std::string>& _nuclides_in_dep_lib);

	// TODO: consider how to construct xs_score
	void InitTallies(std::vector<std::string>& _xs_score);

	void SetNuclideDensities(const int& _mat_id, std::vector<double>& _nuclide_dens);

	void UnpackTallies();

	std::vector<std::string>& nuclides_in_openmc() {
		return nuclides_in_openmc_;
	}

	std::vector<int>& mats_indices() {
		return mats_indices_in_depletion_;
	}

private:
	openmc::Tally* flux_tally_;
	size_t flux_tally_index_;

	openmc::Tally* xs_tally_;
	size_t xs_tally_index_;

	openmc::Tally* heat_tally_;
	size_t heat_tally_index_;

	openmc::Tally* fy_tally_, *weighted_tally_;
	size_t fy_tally_index_, weighted_tally_index_;

	// store the indices of mats which are depletable
	std::vector<int> mats_indices_in_depletion_;

	// nuclides in both openmc and modec
	std::vector<std::string> nuclides_in_openmc_;
};

#endif

