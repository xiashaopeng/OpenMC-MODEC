#include "fission_yields.h"

std::map<int, FYData> FissionYieldsLib::fission_yields_data_;

void FissionYieldsLib::InsertFYData(int _nuclide_zai, std::map<double, std::map<int, double> > _fy_data) {
    auto fy_iter = fission_yields_data_.find(_nuclide_zai);
    if (fy_iter != fission_yields_data_.end()) {
        // fmt::print("");
        return;
    }

    auto duplicate_fy_data = _fy_data;

    std::vector<double> energies_list;
    std::vector<int> fission_products;
    std::vector<std::vector<double> > fission_yields;
    
    for (auto energies_fy_data:_fy_data) {
        auto energy = energies_fy_data.first;
        auto fy_data = energies_fy_data.second;
        
        for (auto fy_data_ele : fy_data) {
            duplicate_fy_data[energy].insert(fy_data_ele);
        }
    }

    for (auto iter = duplicate_fy_data.begin(); iter != duplicate_fy_data.end(); ++iter) {
        
        energies_list.emplace_back(iter->first);

        auto& fy_data = iter->second;
        std::vector<double> fy_value;
        for (auto& ele_value : fy_data) {
            if (iter == duplicate_fy_data.begin())
                fission_products.emplace_back(ele_value.first);
            fy_value.emplace_back(ele_value.second);
        }

        fission_yields.emplace_back(fy_value);
    }

    FYData fydata;
    fydata.fission_energies_ = energies_list;
    fydata.fission_products_ = fission_products;
    fydata.fission_yields_ = fission_yields;

    fission_yields_data_.insert(std::make_pair(_nuclide_zai, fydata));
};

void FissionYieldsLib::EffectiveFYData(int nuclide_zai, double energy, std::vector<int>& effective_nuclide, std::vector<double>& effective_fission_yields) {
    auto fy_data_iter = fission_yields_data_.find(nuclide_zai);
    if (fy_data_iter == fission_yields_data_.end()) {
        return;
    }
    auto fy_data = fy_data_iter->second;

    effective_nuclide = fy_data.fission_products_;

    if (energy <= fy_data.fission_energies_.front()) {
        effective_fission_yields = fy_data.fission_yields_.front();
    }
    else if (energy >= fy_data.fission_energies_.back()) {
        effective_fission_yields = fy_data.fission_yields_.back();
    }
    else {
        auto energy_iter = std::lower_bound(fy_data.fission_energies_.begin(), fy_data.fission_energies_.end(), energy);
        int index = energy_iter - fy_data.fission_energies_.begin();
        if ((*energy_iter) == energy) {
            effective_fission_yields = fy_data.fission_yields_[index];
        }
        else {
            double coeff_upper = ((*energy_iter) - energy)
                                / ((*energy_iter) - (*(energy_iter - 1)));
            
            auto& vec_lower = fy_data.fission_yields_[index - 1];
            auto& vec_upper = fy_data.fission_yields_[index];

            effective_fission_yields.resize(vec_lower.size());
            for (int i = 0; i < effective_fission_yields.size(); ++i) {
                effective_fission_yields[i] = vec_upper[i] - coeff_upper * (vec_upper[i] - vec_lower[i]);
            }
        }
    }
};

//void FissionYieldsLib::EffectiveFYData(const std::string& _nuclide_name, double energy, std::vector<int>& effective_nuclide, std::vector<double>& effective_fission_yields) {
//    int nuclide_zai = NuclideSystem::NuclideZai(_nuclide_name);
//    
//    auto fy_data_iter = fission_yields_data_.find(nuclide_zai);
//    if (fy_data_iter == fission_yields_data_.end()) {
//        return;
//    }
//    auto fy_data = fy_data_iter->second;
//
//    effective_nuclide = fy_data.fission_products_;
//
//    if (energy <= fy_data.fission_energies_.front()) {
//        effective_fission_yields = fy_data.fission_yields_.front();
//    }
//    else if (energy >= fy_data.fission_energies_.back()) {
//        effective_fission_yields = fy_data.fission_yields_.back();
//    }
//    else {
//        auto energy_iter = std::lower_bound(fy_data.fission_energies_.begin(), fy_data.fission_energies_.end(), energy);
//        int index = energy_iter - fy_data.fission_energies_.begin();
//        if ((*energy_iter) == energy) {
//            effective_fission_yields = fy_data.fission_yields_[index];
//        }
//        else {
//            double coeff_upper = ((*energy_iter) - energy)
//                                / ((*energy_iter) - (*(energy_iter - 1)));
//            
//            auto& vec_lower = fy_data.fission_yields_[index - 1];
//            auto& vec_upper = fy_data.fission_yields_[index];
//
//            effective_fission_yields.resize(vec_lower.size());
//            for (int i = 0; i < effective_fission_yields.size(); ++i) {
//                effective_fission_yields[i] = vec_upper[i] - coeff_upper * (vec_upper[i] - vec_lower[i]);
//            }
//        }
//    }
//};
