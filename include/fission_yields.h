#ifndef FISSION_YIELDS_H_
#define FISSION_YIELDS_H_

#include <vector>
#include <string>
#include <map>
// #include "nuclide_system.h"

struct FYData {
    std::vector<double> fission_energies_;
    std::vector<int> fission_products_;
    std::vector<std::vector<double> > fission_yields_;
};

class FissionYieldsLib {
public:
    static void InsertFYData(int _nuclide_zai, std::map<double, std::map<int, double> > _fy_data);

    static void EffectiveFYData(int nuclide_zai, double energy, std::vector<int>& effective_nuclide, std::vector<double>& effective_fission_yields);

    //static void EffectiveFYData(const std::string& _nuclide_name, double energy, std::vector<int>& effective_nuclide, std::vector<double>& effective_fission_yields);

private:
    static std::map<int, FYData> fission_yields_data_;
};

#endif