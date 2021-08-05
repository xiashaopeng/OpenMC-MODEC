/*
 * @Description: class of nuclear data 
 * @version: 
 * @Author: Shaopeng Xia
 * @Date: 2021-07-24 22:56:38
 * @LastEditors: Shaopeng Xia
 * @LastEditTime: 2021-07-27 13:40:37
 */

#ifndef NUCLIDE_DATA_H_
#define NUCLIDE_DATA_H_

#include <vector>
#include <string>

typedef int INDEX;

class XsData {

public:
    double& fission_xs() {return fission_xs_;}
    double& fission_energy() {return fission_energy_;}
    std::vector<double>& xs_values() {return xs_values_;}
    int& num_xs() {return num_xs_;}

    double xs_value(INDEX _index) { return xs_values_[_index];}

    double tot_xs() {
        double tot_xs = fission_xs_;
        for (int i = 0; i < num_xs_; ++i) {
            tot_xs += xs_values_[i];
        }
        return tot_xs;
    };

    // set xs value by the index of the xs type
    void set_xs_value (INDEX _index, const double& _xs_value) {
        if (_index < 0) return;
        if (xs_values_.size() < _index + 1) xs_values_.resize(_index + 1);

        xs_values_[_index] = _xs_value;
    }
    
    // set xs values in a vector
    void set_xs_values (const std::vector<double>& _xs_values) {
        xs_values_ = _xs_values;
    }

    // get xs values for assembling burnup matrix
    void XsValues(std::vector<double>& _xs_values) {
        _xs_values = xs_values_;
        _xs_values.emplace_back(-tot_xs());
    }

    // get xs values for assembling burnup matrix
    void FissionYieldsValues(std::vector<double>& _fission_yields_values) {
        for (int i = 0; i < _fission_yields_values.size(); ++i) {
            _fission_yields_values[i] *= fission_xs_;
        }
    }

private:
    // the value of fission xs
    double fission_xs_{0};

    // neutron energy causing fission
    double fission_energy_{0.0253};

    // the values of other xs
    std::vector<double> xs_values_;

    // the number of total xs
    int num_xs_{0};
};

#endif