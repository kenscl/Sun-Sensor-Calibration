#include "lut.h"
#include <cstdio>

LUT::LUT(std::vector<double> gt_data, std::vector<double> measurement, int num_points) {
    if (gt_data.size() != measurement.size()) {
        throw std::logic_error("datasets must be the same size for LUT");
    }
    std::vector<double> correction_data;
    std::vector<double> meas;
    for (int i = 0; i < gt_data.size(); ++i) {
            correction_data.push_back(gt_data.at(i) - measurement.at(i));
    }
    correction_data = select_evenly_spaced_points(correction_data, num_points);
    meas = select_evenly_spaced_points(measurement, num_points);
    this->parameters = correction_data;
    this->measurement = meas;
    //printf("number of stored points: %zu with %lu bytes\n", this->parameters.size(), this->parameters.size() * sizeof(double));
}

double LUT::at(double x) {
    if (x < this->measurement.at(0)) {
        return x + this->parameters.at(0);
    }
    for (int i = 0; i < this->measurement.size() - 1; ++i) {
        if (x >= this->measurement.at(i) && x < this->measurement.at(i+1)) {
            double d = this->measurement.at(i + 1) - this->measurement.at(i);
            double weigth_left = this->parameters.at(i) * (this->measurement.at(i + 1) - x) / d; 
            double weigth_right = this->parameters.at(i + 1) * ((x - this->measurement.at(i)) / d);
            return x + weigth_left + weigth_right;
        }
    }
    return x + this->parameters.at(this->parameters.size() - 1);
}

std::vector<double> LUT::calc(std::vector<double> x) {
    std::vector<double> ret;
    for (uint i = 0; i < x.size(); i++) {
        ret.push_back(this->at(x.at(i)));
    }
    return ret;
}
