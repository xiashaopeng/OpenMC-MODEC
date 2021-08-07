#include "modec_functions.h"

void MODEC::InitModec(const std::string& _lib_filename) {
    // First initialize decay info
    int nuclide_zai, target_zai, sf_product_zai;
    double nuclide_atomic_mass, nuclide_half_life, nuclide_decay_energy;
    double decay_branching_ratio, sf_yield;
    std::string nuclide_name;
    std::string decay_mode;
    std::string delimiter = ",";

    std::map<int, Nuclide> nuclide_list;

    tinyxml2::XMLDocument lib_file;
    while (lib_file.LoadFile(_lib_filename.c_str())) {
        int error_id = lib_file.ErrorID();
        if (error_id == 0) break;

        if (error_id != 3) {
            std::string info = "ERROR OCCURRED WHEN INPUT FILE WAS BEING READ!\nError MESSAGE: " + std::string(lib_file.ErrorName());
            spdlog::error(info);
            exit(0);
        }
    }

    tinyxml2::XMLElement* modec = lib_file.RootElement();
    tinyxml2::XMLElement* decay_ele = modec->FirstChildElement();

    tinyxml2::XMLElement* nuclide_ele = decay_ele->FirstChildElement();

    while (nuclide_ele) {
        nuclide_zai = atoi(nuclide_ele->Attribute("zai"));
        nuclide_name = nuclide_ele->Attribute("name");
        nuclide_atomic_mass = atof(nuclide_ele->Attribute("atomic_mass"));

        // create a nuclide 
        Nuclide nuclide(nuclide_name, nuclide_zai);
        nuclide.nuclide_atom_mass() = nuclide_atomic_mass;

        // insert new nuclide and check if existed
        auto ret = nuclide_list.insert(std::pair<int, Nuclide>(nuclide_zai, nuclide));
        if (!ret.second) {
            spdlog::warn("The nuclide {} is alreadly existed", nuclide_name);
            nuclide_ele = nuclide_ele->NextSiblingElement();
            continue;
        }

        auto& nuclide_in = nuclide_list[nuclide_zai];

        nuclide_half_life = atof(nuclide_ele->Attribute("half_life"));
        nuclide_decay_energy = atof(nuclide_ele->Attribute("decay_energy"));

        if (nuclide_half_life == 0) {
            nuclide_ele = nuclide_ele->NextSiblingElement();
            continue;
        }

        // set decay lambda and decay heat coefficient
        nuclide_in.set_nuclide_decay_lambda(log(2) / nuclide_half_life);
        nuclide_in.set_nuclide_decay_heat_coeff(nuclide_decay_energy);

        tinyxml2::XMLElement* decay_mode_ele = nuclide_ele->FirstChildElement();
        std::vector<std::string> decay_modes;
        std::vector<int> decay_daughters;
        std::vector<double> decay_branches;

        while (decay_mode_ele) {
            decay_mode = decay_mode_ele->Attribute("mode");
            target_zai = atoi(decay_mode_ele->Attribute("target_zai"));
            decay_branching_ratio = atof(decay_mode_ele->Attribute("branching_ratio"));

            decay_modes.emplace_back(decay_mode);
            if (target_zai == nuclide_zai) {
                // spontaneous fission yields
                tinyxml2::XMLElement* sf_yields = decay_mode_ele->FirstChildElement();
                while (sf_yields) {
                    decay_daughters.emplace_back(atoi(sf_yields->Attribute("products")));
                    decay_branches.emplace_back(decay_branching_ratio * atof(sf_yields->Attribute("yields")));

                    sf_yields = sf_yields->NextSiblingElement();
                }
            }
            else {
                decay_daughters.emplace_back(target_zai);
                decay_branches.emplace_back(decay_branching_ratio);
            }

            decay_mode_ele = decay_mode_ele->NextSiblingElement();
        }

        nuclide_in.set_decay_infos(decay_modes, decay_daughters, decay_branches);

        nuclide_ele = nuclide_ele->NextSiblingElement();
    }

    tinyxml2::XMLElement* xs_ele = decay_ele->NextSiblingElement();
    nuclide_ele = xs_ele->FirstChildElement();

    double reaction_energy;
    std::string reaction_type;

    std::map<int, std::string> fy_subsitutes;

    while (nuclide_ele) {
        nuclide_zai = atoi(nuclide_ele->Attribute("zai"));
        nuclide_name = nuclide_ele->Attribute("name");

        // create a nuclide 
        Nuclide nuclide(nuclide_name, nuclide_zai);

        // insert new nuclide and check if existed
        auto ret = nuclide_list.insert(std::pair<int, Nuclide>(nuclide_zai, nuclide));
        if (!ret.second) {
            spdlog::warn("The nuclide {} is alreadly existed", nuclide_name);
        }

        auto& nuclide_in = nuclide_list[nuclide_zai];

        tinyxml2::XMLElement* nuclide_reactions = nuclide_ele->FirstChildElement();

        std::vector<std::string> xs_types;
        std::vector<int> target_zais;
        std::vector<double> xs_heat_coeffs;
        // branch ratio for 0,1
        std::multimap<std::string, std::pair<int, double> > xs_branch;

        bool has_fy_data{ false };
        std::string fission_nuclide;
        double fission_heat;

        // TODO: 
        while (nuclide_reactions) {
            std::string reac_name = nuclide_reactions->Name();
            if (reac_name == "reaction") {
                if (nuclide_reactions->Attribute("target_zai")) {
                    target_zai = atoi(nuclide_reactions->Attribute("target_zai"));
                    target_zais.emplace_back(target_zai);

                    reaction_type = nuclide_reactions->Attribute("type");
                    xs_types.emplace_back(reaction_type);

                    xs_heat_coeffs.emplace_back(atof(nuclide_reactions->Attribute("Q")) * 1.0e-6); // convert eV to MeV

                    if (nuclide_reactions->Attribute("branching_ratio")) {
                        double branching_ratio_xs = atof(nuclide_reactions->Attribute("branching_ratio"));
                        xs_branch.insert(std::make_pair(reaction_type, std::make_pair(target_zai, branching_ratio_xs)));
                    }
                }
                else {
                    nuclide_in.set_fission_heat(atof(nuclide_reactions->Attribute("Q")) * 1.0e-6); // W -> MW
                }
            }
            else {
                if (nuclide_reactions->Attribute("substitute")) {
                    fy_subsitutes.insert(std::make_pair(nuclide_zai, nuclide_reactions->Attribute("substitute")));
                    //nuclide_in.fission_nuclide() = nuclide_reactions->Attribute("substitute");
                }
                else {
                    // FIXME: should use nuclide name not nuclide zai
                    nuclide_in.fission_nuclide() = nuclide_zai;

                    tinyxml2::XMLElement* neutron_fission_energy_ele = nuclide_reactions->FirstChildElement();
                    std::map<double, std::map<int, double> > nfy_map;
                    while (neutron_fission_energy_ele) {
                        double fission_neutron_energy = atof(neutron_fission_energy_ele->Attribute("energy"));
                        tinyxml2::XMLElement* nfy_ele = neutron_fission_energy_ele->FirstChildElement();

                        std::map<int, double> nf_yield_map;
                        while (nfy_ele) {
                            int product = atoi(nfy_ele->Attribute("product"));
                            double yield = atof(nfy_ele->Attribute("yields"));

                            nf_yield_map.insert(std::map<int, double>::value_type(product, yield));

                            nfy_ele = nfy_ele->NextSiblingElement();
                        }
                        nfy_map.insert(std::map<double, std::map<int, double> >::value_type(fission_neutron_energy, nf_yield_map));

                        neutron_fission_energy_ele = neutron_fission_energy_ele->NextSiblingElement();
                    }

                    nuclide_in.InsertFissionYieldsLib(nfy_map);
                }
            }

            nuclide_in.set_xs_infos_except_fission(xs_types, target_zais, xs_heat_coeffs, xs_branch);
            nuclide_reactions = nuclide_reactions->NextSiblingElement();
        }

        nuclide_ele = nuclide_ele->NextSiblingElement();
    }

    NuclideSystem::NuclideHelper(nuclide_list, fy_subsitutes);
};

void MODEC::InitDepZones(DepZones& _dep_zones, std::vector<int>& _zone_indices) {
    //std::vector<int32_t> material_indices;
    //std::vector<int> a;

    //// iterate over mats
    //for (auto& entry : openmc::model::material_map) {
    //    material_indices.push_back(entry.first);

    //    int* nuclide_indices;
    //    double* nuclide_dens;
    //    const char* nuclide_names;
    //    int num_nucl;

    //    // get initial densities using material id 
    //    openmc_material_get_densities(entry.first, &nuclide_indices, &nuclide_dens, &num_nucl);

    //    for (int i = 0; i < num_nucl; ++i) {
    //        openmc_nuclide_name(i, &nuclide_names);
    //    }
    //    
    //}
    // _dep_zones.Initialize(material_indices);
    _dep_zones.Initialize(_zone_indices);
}

SparseMatrix<double> MODEC::AssembleMat(NuclideSystem& _nuclide_system, const double& norm_factor) {
    int num_nucl = Nuclide::TotalNuclideNum();

    // construct decay matrix
    SparseMatrix<double> decay_mat(num_nucl);
    _nuclide_system.AssembleDecayMat(decay_mat.row_indices_, decay_mat.col_compressed_slice_, decay_mat.values_);
    decay_mat.DuplicateElements();

    // construct xs matrix
    SparseMatrix<double> xs_mat(num_nucl);
    _nuclide_system.AssembleDecayMat(xs_mat.row_indices_, xs_mat.col_compressed_slice_, xs_mat.values_);
    xs_mat.DuplicateElements();

    return decay_mat + xs_mat * norm_factor;
};

void MODEC::CalcPhiFunction(SparseMatrix<double> _burnup_matrix, SymbolicSpMat& sym_spmat, std::map<int, std::vector<double> >& _vec_b) {

}