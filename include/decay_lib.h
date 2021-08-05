#ifndef DECAY_LIB_H_
#define DECAY_LIB_H_

#include <vector>
#include <string>

class DecayLib {
public:
    double& decay_lambda() {return decay_lambda_;}
    double& decay_heat_coeff() {return decay_heat_coeff_;}
    double& decay_radiotoxicity_coeff() {return decay_radiotoxicity_coeff_;}
    double& decay_ampcs_coeff() {return decay_ampcs_coeff_;}
    double& decay_wmpcs_coeff() {return decay_wmpcs_coeff_;}

    std::vector<std::string>& decay_mode() {return decay_modes_;}
    std::vector<int>& decay_daughters() {return decay_daughters_;}
    std::vector<double>& decay_branches() {return decay_branches_;}

private:
    double decay_lambda_{0};
    double decay_heat_coeff_{0};
    double decay_radiotoxicity_coeff_{0};
    double decay_ampcs_coeff_{0};
    double decay_wmpcs_coeff_{0};
    //std::map<std::string, std::map<int, double> > decay_branch_; // key is decay mode; value map is the pair of daughter and the branch value
    std::vector<std::string> decay_modes_;
    std::vector<int> decay_daughters_;
    std::vector<double> decay_branches_;
};

#endif