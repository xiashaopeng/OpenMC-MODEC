#ifndef OPENMC_OPERATOR_H_
#define OPENMC_OPERATOR_H_

#include <openmc/simulation.h>

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/material.h"
#include "openmc/message_passing.h"
#include "openmc/summary.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_cell.h"
#include "openmc/tallies/tally.h"

#include "transport_operator.h"
#include "dep_zones.h"

class OpenMCOperator : public TransportOperator {
public:

	void CreateTallies() {
		// create a new cell filter
		auto cell_filter = openmc::Filter::create<openmc::CellFilter>();

		// add all cells to the cell filter
		std::vector<int32_t> cell_indices;
		for (auto& entry : openmc::model::cell_map) {
			cell_indices.push_back(entry.second);
		}
		// sort to make sure the cell bins appear in the same
		// order as the test relying on the openmc exe
		std::sort(cell_indices.begin(), cell_indices.end());
		cell_filter->set_cells(cell_indices);

		// create a flux tally
		flux_tally_ = openmc::Tally::create();
		std::vector<openmc::Filter*> filters = { cell_filter };
		flux_tally_->set_filters(filters);
		flux_tally_->set_scores({ "flux" });

		// add all materials' id to create NuclideSystem
		std::vector<int32_t> material_indices;
		for (auto& entry : openmc::model::material_map) {
			material_indices.push_back(entry.first);
		}

		//openmc_get_nuclide_index


	}

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

private:
	openmc::Tally* flux_tally_;
	openmc::Tally* xs_tally_;
	openmc::Tally* heat_tally_;
	openmc::Tally* fy_tally_;

	// nuclides in both openmc and modec
	std::vector<std::string> nuclides_in_openmc_;
};

#endif

