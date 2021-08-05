#ifndef XS_LIB_H_
#define XS_LIB_H_

#include <vector>
#include <string>
#include <map>

typedef int INDEX;

class XsLib {
public:
    // expose all member variables to Nuclide class
    std::map<std::string, INDEX>& xs_types() {return xs_types_;}
    std::map<std::string, INDEX>& xs_bytypes() {return xs_bytypes_;}
    std::multimap<std::string, std::pair<int, double> >& xs_branches() {return xs_branches_;}

    std::vector<int>& xs_products() {return xs_products_;}
    std::vector<double>& xs_heat_coeffs() {return xs_heat_coeffs_;}

    double& fission_heat_coeffs() {return fission_heat_coeffs_;}

    // functional methods

    // return index of a type of xs
    void get_xs_index (const std::string& _xs_type, INDEX& _prod_index, INDEX& _by_prod_index) {
        auto iter = xs_types_.find(_xs_type);
        if (iter == xs_types_.end()) {
            _prod_index = -1;
        }
        else {
            _prod_index = iter->second;
        }
        
        iter = xs_bytypes_.find(_xs_type);
        if (iter == xs_types_.end()) {
            _by_prod_index = -1;
        }
        else {
            _by_prod_index = iter->second;
        }
    }

    // return heat coefficient
    void get_xs_heat_coeff (const std::string& _xs_type, double& heat_coeffs) {
        auto iter = xs_types_.find(_xs_type);
        if (iter == xs_types_.end()) {
            heat_coeffs = 0;
            return;
        }

        heat_coeffs = xs_heat_coeffs_[iter->second];
    }

private:
    // store all the reactions of the nuclide except the fission reaction
    std::map<std::string, INDEX> xs_types_;
    // store the by-products of all the reactions
    std::map<std::string, INDEX> xs_bytypes_;
    
    // store branch ratio
    std::multimap<std::string, std::pair<int, double> > xs_branches_;

    std::vector<int> xs_products_;
    std::vector<double> xs_heat_coeffs_;

    // store fission heat coefficient seperately
    double fission_heat_coeffs_{0};
};

#endif