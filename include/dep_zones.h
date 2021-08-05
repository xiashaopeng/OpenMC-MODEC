#ifndef DEP_ZONES_H_
#define DEP_ZONES_H_

#include "nuclide_system.h"
#include "sparse_lu.h"

// xs data buffer for data transfer
struct XsDataBuffer {
    XsDataBuffer(int& _nuclide_zai, std::vector<std::string>& _xs_types, std::vector<double>& _xs_value)
        : nuclide_zai_(_nuclide_zai), xs_types_(_xs_types), xs_value_(_xs_value) {};

    int& nuclide_zai_;
    std::vector<std::string>& xs_types_;
    std::vector<double>& xs_value_;
};

// multi depletion zones
// TODO: calc_helper is considered remove from this class
class DepZones {
public:
    // singleton pattern to construct dep zones
    static DepZones& CreateDepZones();

    // initialize dep_zones instance
    void Initialize(const std::vector<int32_t>&);

    std::vector<double>& nuclide_densities(const int& _zone_id) {
        auto& nuclide_system = burnup_zones_[_zone_id];
        return nuclide_system.nuclide_densities();
    }

    double& power_dens(const int& _zone_id) {
        return power_dens_[_zone_id];
    }

    double& neutron_flux(const int& _zone_id) {
        return neutron_fluxes_[_zone_id];
    }

    double normalization_factor(const int& _zone_id) {
        if (modes_[_zone_id] == "flux") {
            return neutron_fluxes_[_zone_id];
        }
        else if (modes_[_zone_id] == "power") {
            return power_dens_[_zone_id];
        }
    }

    std::vector<double>& time_steps() {
        return time_steps_;
    }

    // get burnup matrix for a depletion zone
    SparseMatrix<double> burnup_matrix(const int& _zone_id);

    // update xs values for a single depletion zone
    void UpdateXsPerZone(const int& _zone_id, std::vector<XsDataBuffer>& _xs_data_buffer);

    // update xs values
    void UpdateXs(std::unordered_map<int, std::vector<XsDataBuffer> >& _xs_data_buffers);

private:
    DepZones() {};
    ~DepZones() {};
    DepZones(const DepZones&);
    DepZones& operator=(const DepZones&);

    // key is the id of the zone
    std::unordered_map<int, NuclideSystem> burnup_zones_;

    // power densities of each burnup zone
    std::unordered_map<int, double> power_dens_;

    // neutron flux of each burnup zone
    std::unordered_map<int, double> neutron_fluxes_;

    // calculation mode : constant power, constant flux
    std::unordered_map<int, std::string> modes_;

    // burnup time steps for simulations
    std::vector<double> time_steps_;

    // symbolic lu factorization helper
    std::unordered_map<int, SymbolicSpMat> calc_helper;

    /* void Simulate(int _zone_id) {
        auto iter_nuclide_system = burnup_zones_.find(_zone_id);
        if (iter_nuclide_system == burnup_zones_.end()) {
            return;
        }
        auto& nuclide_system = iter_nuclide_system->second;

        // set initial densities
        auto& init_dens = nuclide_system.nuclide_densities();

        // get normalization factor for burnup matrix
        double norm_factor;
        if (modes_[_zone_id] == "flux") {
            norm_factor = neutron_fluxes_[_zone_id];
        }
        else if (modes_[_zone_id] == "power") {
            norm_factor = power_dens_[_zone_id];
        }

        // get burnup matrix of the depletion zone
        SparseMatrix<double> mat(MODEC::AssembleMat(nuclide_system, norm_factor));

        // simulate burnup equations
        auto& sp_mat = calc_helper[_zone_id];

    }; */


};


#endif

