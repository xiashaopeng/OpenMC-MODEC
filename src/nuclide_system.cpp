#include "nuclide_system.h"

std::unordered_map<INDEX, Nuclide> NuclideSystem::nuclides_helper_;
std::unordered_map<int, INDEX> NuclideSystem::zai_helper_;
MODEC::bimap<int, std::string> NuclideSystem::name_helper_; 

void NuclideSystem::NuclideHelper(const std::map<int, Nuclide>& _nuclides,
    const std::map<int, std::string>& _fy_subsitute) {

    // set nuclides helper
    int cnt = 0;
    for (auto& ele:_nuclides) {
        nuclides_helper_.insert(std::make_pair(cnt, ele.second));
        zai_helper_.insert(std::make_pair(ele.first, cnt));
        name_helper_.insert(std::make_pair(ele.first, ele.second.nuclide_name()));
    }

    // update the fission nuclide infos
    for (auto& fy_ele : _fy_subsitute) {
        int target_zai = fy_ele.first;
        INDEX index = zai_helper_[target_zai];

        std::string _name = fy_ele.second;
        int zai = name_helper_.get_key(_name);

        nuclides_helper_[index].fission_nuclide() = zai;
    }
};

int NuclideSystem::NuclideZai(const std::string& _name) {
    return name_helper_.get_key(_name);
};

std::string NuclideSystem::NuclideName(const int& _zai) {
    return name_helper_.get_value(_zai);
};

void NuclideSystem::InitXsData() {
    int num_nucl = Nuclide::TotalNuclideNum();
    xs_data_.resize(num_nucl);

    for (auto& nuclide : nuclides_helper_) {
        xs_data_[nuclide.first].num_xs() = nuclide.second.num_xs();
    }
};

void NuclideSystem::AssembleDecayMat(std::vector<int>& _row_indices, std::vector<int>& _col_compressed_slice, std::vector<double>& _value) {
    int num_nucl = Nuclide::TotalNuclideNum();

    std::vector<int> row_indices;
    std::vector<int> col_compressed_slice;

    std::vector<double> value;

    int cnt_element = 0;
    for (int i = 0; i < num_nucl; ++i) {
        auto& nuclide = nuclides_helper_[i];

        std::vector<int> target_zai;
        std::vector<double> decay_values;
        nuclide.GetDecayValues(target_zai, decay_values);
        
        col_compressed_slice.emplace_back(cnt_element);
        for (int j = 0; j < target_zai.size(); ++j) {
            row_indices.emplace_back(zai_helper_[target_zai[j]]);
            value.emplace_back(decay_values[j]);
        }

        cnt_element += decay_values.size();

        // add remove coefficient
        if (remove_coeffs_.size() > 0 && remove_coeffs_[i] != 0) {
            decay_values.back() -= remove_coeffs_[i];
        }
    }
    col_compressed_slice.emplace_back(cnt_element);

}

void NuclideSystem::AssembleXsMat(std::vector<int>& _row_indices, std::vector<int>& _col_compressed_slice, std::vector<double>& _value) {
    int num_nucl = Nuclide::TotalNuclideNum();

    std::vector<int> row_indices;
    std::vector<int> col_compressed_slice;

    std::vector<double> value;

    int cnt = 0;
    int cnt_element = 0;
    for (int i = 0; i < num_nucl; ++i) {
        auto& nuclide = nuclides_helper_[i];

        // construct matrix using non-fission xs
        std::vector<int> target_zai;
        std::vector<double> xs_values;
        nuclide.GetXsDaughters(target_zai);
        xs_data_[i].XsValues(xs_values);
        
        col_compressed_slice.emplace_back(cnt_element);
        for (int j = 0; j < target_zai.size(); ++j) {
            row_indices.emplace_back(zai_helper_[target_zai[j]]);
            value.emplace_back(xs_values[j]);
        }

        // construct matrix using fission xs
        std::vector<int> fission_products;
        std::vector<double> fission_values;
        double energy = xs_data_[i].fission_energy();

        nuclide.CalcEffectiveFissionYields(energy, fission_products, fission_values);
        xs_data_[i].FissionYieldsValues(fission_values);
        for (int j = 0; j < fission_products.size(); ++j) {
            row_indices.emplace_back(zai_helper_[fission_products[j]]);
            value.emplace_back(fission_values[j]);
        }

        cnt_element += xs_values.size();
        cnt_element += fission_products.size();
    }
    col_compressed_slice.emplace_back(cnt_element);
}

void NuclideSystem::UpdateXsValues(const int& _nuclide_zai, const std::vector<std::string>& _xs_types, const std::vector<double>& _xs_value) {
    INDEX index = zai_helper_[_nuclide_zai];
    auto& nuclide = nuclides_helper_[index];
    auto& xs_data = xs_data_[index];

    for (int i = 0; i < _xs_types.size(); ++i) {
        if (_xs_types[i] == "fission") {
            xs_data.fission_xs() = _xs_value[i];
            continue;
        }

        INDEX prod_index{-1}, by_prod_index{-1};

        nuclide.GetXsIndex(_xs_types[i], prod_index, by_prod_index);
        xs_data.set_xs_value(prod_index, _xs_value[i]);
        xs_data.set_xs_value(by_prod_index, _xs_value[i]);

        // set values for the xs which has branches
        std::vector<INDEX> indices;
        std::vector<double> branch_ratios;
        nuclide.GetBranchesIndices(_xs_types[i], indices, branch_ratios);
        for (int j = 0; j < indices.size(); ++j) {
            xs_data.set_xs_value(indices[j], _xs_value[i] * branch_ratios[j]);
        }
    }
} 

void NuclideSystem::CalcPowerFluxKernel() {
    int num_nucl = Nuclide::TotalNuclideNum();

    kernel_ = 0;

    for (int i = 0; i < num_nucl; ++i) {
        double& fission_xs = xs_data_[i].fission_xs();
        if (fission_xs == 0) continue;

        double& dens = nuclide_densities_[i];
        double& fission_heat_coeffs = nuclides_helper_[i].fission_heat();

        kernel_ += fission_xs * dens * fission_heat_coeffs;
    }

    // calculate kernel using fission xs and capture xs
    for (int i = 0; i < num_nucl; ++i) {
        double& fission_xs = xs_data_[i].fission_xs();
        double& fission_heat_coeffs = nuclides_helper_[i].fission_heat();

        int cap_index{-1}, no_use_index{-1};
        nuclides_helper_[i].GetXsIndex("(n,gamma)", cap_index, no_use_index);
        double gamma_xs = xs_data_[i].xs_value(cap_index);

        double gamma_heat_coeffs;
        nuclides_helper_[i].GetXsHeatCoeff("(n,gamma)", gamma_heat_coeffs);

        double& dens = nuclide_densities_[i];

        kernel_ += (fission_xs * fission_heat_coeffs + gamma_xs * gamma_heat_coeffs) * dens;
    }

    // calculate kernel using all heat coeffs
    for (int i = 0; i < num_nucl; ++i) {
        double& fission_xs = xs_data_[i].fission_xs();
        double& fission_heat_coeffs = nuclides_helper_[i].fission_heat();

        int cap_index{-1}, no_use_index{-1};
        nuclides_helper_[i].GetXsIndex("(n,gamma)", cap_index, no_use_index);
        double gamma_xs = xs_data_[i].xs_value(cap_index);

        double gamma_heat_coeffs;
        nuclides_helper_[i].GetXsHeatCoeff("(n,gamma)", gamma_heat_coeffs);

        double& dens = nuclide_densities_[i];

        kernel_ += (fission_xs * fission_heat_coeffs + gamma_xs * gamma_heat_coeffs) * dens;
    }

};
