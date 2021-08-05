#include "dep_zones.h"

DepZones& DepZones::CreateDepZones() {
    static DepZones instance;
    return instance;
};

void DepZones::Initialize(const std::vector<int32_t>& _zone_incices) {
    for (int i = 0; i < _zone_incices.size(); ++i) {
        burnup_zones_.insert(std::make_pair(i, NuclideSystem()));
        power_dens_.insert(std::make_pair(i, 0.0));
        neutron_fluxes_.insert(std::make_pair(i, 0.0));
        modes_.insert(std::make_pair(i, "power"));
        calc_helper.insert(std::make_pair(i, SymbolicSpMat()));
    }
};

SparseMatrix<double> DepZones::burnup_matrix(const int& _zone_id) {
    auto& nuclide_system = burnup_zones_[_zone_id];

    int num_size = Nuclide::TotalNuclideNum();

    // construct decay matrix
    SparseMatrix<double> decay_mat(num_size);
    nuclide_system.AssembleDecayMat(decay_mat.row_indices_, decay_mat.col_compressed_slice_, decay_mat.values_);
    decay_mat.DuplicateElements();

    // construct xs matrix
    SparseMatrix<double> xs_mat(num_size);
    nuclide_system.AssembleXsMat(xs_mat.row_indices_, xs_mat.col_compressed_slice_, xs_mat.values_);
    xs_mat.DuplicateElements();

    // get normalization factor
    auto norm_factor = normalization_factor(_zone_id);

    return decay_mat + xs_mat * norm_factor;
}

void DepZones::UpdateXsPerZone(const int& _zone_id, std::vector<XsDataBuffer>& _xs_data_buffer) {
    for (auto& ele : _xs_data_buffer) {
        burnup_zones_[_zone_id].UpdateXsValues(ele.nuclide_zai_, ele.xs_types_, ele.xs_value_);
    }
}

void DepZones::UpdateXs(std::unordered_map<int, std::vector<XsDataBuffer> >& _xs_data_buffers) {
    for (auto& ele : _xs_data_buffers) {
        UpdateXsPerZone(ele.first, ele.second);
    }
}