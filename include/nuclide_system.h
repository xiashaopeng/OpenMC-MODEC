#ifndef NUCLIDE_SYSTEMS_H_
#define NUCLIDE_SYSTEMS_H_

#include <unordered_map>

#include "nuclide.h"
#include "xs_data.h"
#include "bimap.h"

//TODO: should consider the polynominal approximation to xs
// store nuclide's info
class NuclideSystem {
public:
    // initialize nuclide system in constructor
    NuclideSystem() {
        this->InitXsData();
    };

    // static method for initializing helpers
    static void NuclideHelper(const std::map<int, Nuclide>& _nuclides,
        const std::map<int, std::string>& _fy_subsitute);

    // get nuclide zai using its name
    static int NuclideZai(const std::string& _name);

    // get nuclide name using its zai
    static std::string NuclideName(const int& _zai);

    // initialize xs data using nuclide helper
    void InitXsData();

    // expose nuclide densities for using and assigning
    std::vector<double>& nuclide_densities() {return nuclide_densities_;};

    // expose remove coeffs
    std::vector<double>& remove_coeffs() {return remove_coeffs_;};

    // update xs values for a single nuclide
    void UpdateXsValues(const int& _nuclide_zai, const std::vector<std::string>& _xs_types, const std::vector<double>& _xs_value);

    //void AssembleMat(std::vector<int>& _row_indices, std::vector<int>& _col_compressed_slice, std::vector<double>& _value);

    // construct decay matrix
    // TODO:should consider add remove coefficients
    void AssembleDecayMat(std::vector<int>& _row_indices, std::vector<int>& _col_compressed_slice, std::vector<double>& _value);

    // construct xs matrix
    void AssembleXsMat(std::vector<int>& _row_indices, std::vector<int>& _col_compressed_slice, std::vector<double>& _value);

private:
    // nuclide densities used for burnup calculation
    std::vector<double> nuclide_densities_;

    // nuclide xs data used for burnup calculation
    std::vector<XsData> xs_data_;

    // remove coefficients for each nuclide
    // zero size if no removable nuclides
    std::vector<double> remove_coeffs_;

    // flux and power for the nuclide system
    double kernel_;

    // the key is nuclide's index in the vector of densities
    static std::unordered_map<INDEX, Nuclide> nuclides_helper_;
    
    // the key is nuclide's zai; 
    // the value is the index of nuclide in the vector of densities
    static std::unordered_map<int, INDEX> zai_helper_;

    // bimap to map zai and name 
    static MODEC::bimap<int, std::string> name_helper_; 

    void CalcPowerFluxKernel();
};


#endif