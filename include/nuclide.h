#ifndef NUCLIDES_H_
#define NUCLIDEs_H_

#include <string>
#include <map>
#include <vector>
#include <sstream>

#include "decay_lib.h"
#include "xs_lib.h"
#include "fission_yields.h"

template<typename T>
class NuclidesNumCounter {
public:
    // static int num_nuclides() {
    //     return objects_alive;
    // }

protected:
    NuclidesNumCounter()
    {
        ++objects_created;
        ++objects_alive;
    }
    
    NuclidesNumCounter(const NuclidesNumCounter&)
    {
        ++objects_created;
        ++objects_alive;
    }

    ~NuclidesNumCounter() // objects should never be removed through pointers of this type
    {
        --objects_alive;
    }    

    static int objects_created;
    static int objects_alive;
};

template <typename T> int NuclidesNumCounter<T>::objects_created( 0 );
template <typename T> int NuclidesNumCounter<T>::objects_alive( 0 );

class Nuclide : public NuclidesNumCounter<Nuclide> {
public:
    // constructors
    Nuclide() {};
    Nuclide(const std::string& _nuclide_name): name_(_nuclide_name) {};
    Nuclide(int _zai): zai_(_zai) {};
    Nuclide(const std::string& _nuclide_name, int _zai): name_(_nuclide_name), zai_(_zai) {};

    bool operator<(const Nuclide& rhs) const {
        return this->zai_ < rhs.nuclide_zai();
    }

    bool operator==(const Nuclide& rhs) const {
        return this->zai_ == rhs.nuclide_zai();
    }

    // return the number of all nuclides 
    static int TotalNuclideNum() {
        return objects_alive;
    }

    // set nuclide's name
    void set_nuclide_name(const std::string& _nuclide_name) const { name_ = _nuclide_name; }
    // set nuclide's zai
    void set_nuclide_zai(int _nuclide_zai) const { zai_ = _nuclide_zai; }
    // set nuclide's atom mass
    void set_nuclide_atom_mass(double _atom_mass) const { atom_mass_ = _atom_mass; }


    // set nuclide's decay lambda
    void set_nuclide_decay_lambda(double _decay_lambda) const { decay_lib_.decay_lambda() = _decay_lambda; }
    // set nuclide's decay heat coefficient
    void set_nuclide_decay_heat_coeff(double _decay_heat_coeff) const { decay_lib_.decay_heat_coeff() = _decay_heat_coeff; }
    // set nuclide's decay radiotoxicity coefficient
    void set_nuclide_decay_radiotoxicity_coeff(double _decay_radiotoxicity_coeff) const { decay_lib_.decay_radiotoxicity_coeff() = _decay_radiotoxicity_coeff; }
    // set nuclide's decay ampc coefficient
    void set_nuclide_decay_ampc_coeff(double _decay_ampcs_coeff) const { decay_lib_.decay_ampcs_coeff() = _decay_ampcs_coeff; }
    // set nuclide's decay wmpc coefficient
    void set_nuclide_decay_wmpc_coeff(double _decay_wmpcs_coeff) const { decay_lib_.decay_wmpcs_coeff() = _decay_wmpcs_coeff; }
    // set nuclide's decay infos including decay modes, decay daughters and decay branches.
    void set_decay_infos(std::vector<std::string>& _decay_modes, std::vector<int>& _decay_daughters, std::vector<double>& _decay_branches) const {
        // find by-products
        for (int i = 0; i < _decay_modes.size(); ++i) {
            std::istringstream ss_decay_mode(_decay_modes[i]);
            std::string buffer;
            while(getline(ss_decay_mode, buffer, ',')) {
                if (buffer == "alpha") {
                    _decay_daughters.emplace_back(20040);
                    _decay_branches.emplace_back(_decay_branches[i]);
                }
                else if (buffer == "p") {
                    _decay_daughters.emplace_back(10010);
                    _decay_branches.emplace_back(_decay_branches[i]);
                }
            }
        }

        decay_lib_.decay_mode() = _decay_modes;
        decay_lib_.decay_daughters() = _decay_daughters;
        decay_lib_.decay_branches() = _decay_branches;
    }

    void GetDecayValues(std::vector<int>& _target_zai, std::vector<double>& _decay_values) const {
        auto& lambda = decay_lib_.decay_lambda();
        auto decay_target = decay_lib_.decay_daughters();
        auto decay_values = decay_lib_.decay_branches();

        for (int i = 0; i < decay_target.size(); ++i) {
            decay_values[i] *= lambda;
        }

        decay_target.emplace_back(zai_);
        decay_values.emplace_back(-lambda);

        _target_zai = decay_target;
        _decay_values = decay_values;
    }


    // set nuclide's cross section data
    // the orders of all vectors must be the same!!
    void set_xs_infos_except_fission(std::vector<std::string>& _xs_types, std::vector<int>& _target_zais, std::vector<double>& _xs_heat_coeffs, std::multimap<std::string, std::pair<int,double> >& _xs_branches) const {
        
        // the values of the map are the indexes of xs types.
        std::map<std::string, INDEX> xs_types;
        std::map<std::string, INDEX> by_xs_types;
        for (int i = 0; i < _xs_types.size(); ++i) {
            auto ret = xs_types.insert(std::make_pair(_xs_types[i], i));
            if (!ret.second) {
                auto xs_branch = _xs_branches.equal_range(_xs_types[i]);
                for (auto iter = xs_branch.first; iter != xs_branch.second; ++iter) {
                    auto& zai = iter->second;
                    if (zai.first == _target_zais[i]) {
                        zai.first = i;
                    }
                    else {
                        zai.first = xs_types[_xs_types[i]];
                    }
                }

                xs_types[_xs_types[i]] = -1;

                continue;
            }

            // find the alpha or 
            auto sub_str = _xs_types[i].substr(1, _xs_types[i].size() - 2);
            std::istringstream ss_xs_type(sub_str);
            std::string xs_buffer;

            while(getline(ss_xs_type, xs_buffer, ',')) {
                if (xs_buffer == "a") {
                    _target_zais.emplace_back(20040);
                    by_xs_types.insert(std::make_pair(_xs_types[i], _target_zais.size() - 1));
                }
                else if (xs_buffer == "p") {
                    _target_zais.emplace_back(10010);
                    by_xs_types.insert(std::make_pair(_xs_types[i], _target_zais.size() - 1));
                }
                else if (xs_buffer == "d") {
                    _target_zais.emplace_back(10020);
                    by_xs_types.insert(std::make_pair(_xs_types[i], _target_zais.size() - 1));
                }
                else if (xs_buffer == "t") {
                    _target_zais.emplace_back(10030);
                    by_xs_types.insert(std::make_pair(_xs_types[i], _target_zais.size() - 1));
                }
            }
        }

        xs_lib_.xs_types() = xs_types;
        xs_lib_.xs_bytypes() = by_xs_types;

        xs_lib_.xs_products() = _target_zais;
        xs_lib_.xs_heat_coeffs() = _xs_heat_coeffs;

        xs_lib_.xs_branches() = _xs_branches;
    }
    
    // get index of cross secion by the xs type
    void GetXsIndex(const std::string& _xs_type, INDEX& _prod_index, INDEX& _by_prod_index) const {
        xs_lib_.get_xs_index(_xs_type, _prod_index, _by_prod_index);
    }

    // get index of xs with branch ratios
    void GetBranchesIndices(const std::string& _xs_type, std::vector<int>& _indices, std::vector<double>& _branch_ratios) const {
        auto& helper = xs_lib_.xs_branches();
        auto ret = helper.equal_range(_xs_type);
        for (auto iter = ret.first; iter != ret.second; ++iter) {
            _indices.emplace_back(iter->second.first);
            _branch_ratios.emplace_back(iter->second.second);
        }
    }

    // return xs daughters to construct burnup matrix
    void GetXsDaughters(std::vector<int>& _target_zai) const {
        auto target_zai = xs_lib_.xs_products();
        target_zai.emplace_back(zai_);

        _target_zai = target_zai;
    }

    void GetXsHeatCoeff(const std::string& _xs_type, double& _heat_coeff) {
        INDEX index = xs_lib_.xs_types()[_xs_type];
        _heat_coeff = xs_lib_.xs_heat_coeffs()[index];
    }


    // assign nuclide which stores fission yields
    void set_fission_heat(const double& _heat_coeffs) const {
        xs_lib_.fission_heat_coeffs() = _heat_coeffs;
    }

    // set data of fission yields
    void InsertFissionYieldsLib(std::map<double, std::map<int, double> > _fy_data) const {
        FissionYieldsLib::InsertFYData(zai_, _fy_data);
    }

    // return effective fission yields by giving average energy of neutrons 
    // causing fission reactions
    void CalcEffectiveFissionYields(double _energy, std::vector<int>& effective_nuclide, std::vector<double>& effective_fission_yields) const {
        FissionYieldsLib::EffectiveFYData(fission_nuclide_, _energy, effective_nuclide, effective_fission_yields);
    }

    // get nuclide's name
    std::string& nuclide_name() const {return name_;}
    // get nuclide's zai
    int& nuclide_zai() const {return zai_;}
    // get nuclide's atom mass
    double& nuclide_atom_mass() const {return atom_mass_;}    
    // get nuclide's decay lambda
    double& decay_lambda() const {return decay_lib_.decay_lambda();}
    // get the number of nuclide reactions
    int num_xs() const {return xs_lib_.xs_types().size();}
    // get fission heat coeffs
    double& fission_heat() {return xs_lib_.fission_heat_coeffs();}
    // get fission nuclide
    int& fission_nuclide() {return fission_nuclide_;}

    
private:
    mutable std::string name_;
    mutable int zai_;
    mutable double atom_mass_;

    mutable DecayLib decay_lib_;
    mutable XsLib xs_lib_;
    mutable int fission_nuclide_;

    // store nuclide densities at all steps
    // mutable std::vector<double> density_history_;

};

#endif