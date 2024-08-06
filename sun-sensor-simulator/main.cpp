#include <cstdio>
#include <vector>
#include "math/matlib.h"
#include "sun_sensor.h"
#include "sun_sensor_data_parser.h"

// Plotting includes
#include "polynomial_fit/polynomial_fit.h"
#include "sin_fit/sin_fit.h"
#include "lut/lut.h"
#include "matplotlibcpp.h"
#include "sin_fit/extended_sin_fit.h"
#include "sin_fit/extended_sin_fit_alt.h"
#include "lut/c-spline.h"
#include "lut/b-spline.h"

//testing 
#include "test/test.h"

namespace plt = matplotlibcpp;


std::vector<double> filter(std::vector<double> data, int num_points){
    std::vector<double> ret;
    for (int i = 0; i < data.size(); ++i) {
        if (i % (int) (data.size() / (num_points - 1)) == 0) {
            ret.push_back(data.at(i));
        }
    }
    return ret;
}

void filter_45(std::vector<double> *gt, std::vector<double> *meas) {
    std::vector<double> gt_ret, meas_ret;
    for (int i = 0; i < gt->size(); i++) {
        if (gt->at(i) > -45 * D2R && gt->at(i) < 45 * D2R){
            gt_ret.push_back(gt->at(i));
            meas_ret.push_back(meas->at(i));
        }
    }
    *gt = gt_ret;
    *meas = meas_ret;
}

int main () {

    std::vector<std::vector<double>> data, data_valid;
    data = parse_csv("../../Sunsensor_Data/28.08.2021 17.20.54 ADPD2140_X.txt");
    std::vector<double> stage_angle = data.at(0);
    std::vector<double> X = data.at(1);
    std::vector<double> Y = data.at(2);

    std::vector<double> pf_cor;
    std::vector<double> sf_cor;
    std::vector<double> ef_cor;
    std::vector<double> ef_alt_cor;
    std::vector<double> lut_cor;

    data_valid = parse_csv("../../Sunsensor_Data/28.08.2021 16.37.25 ADPD2140_Y.txt");
    std::vector<double> stage_angle_y = data_valid.at(0);
    std::vector<double> Y_y= data_valid.at(1);
    std::vector<double> X_y= data_valid.at(2);

    std::vector<double> error_x, error_y;
    for(int i = 0; i < stage_angle.size(); i++) {
        error_x.push_back(stage_angle.at(i) - X.at(i));
        error_y.push_back(stage_angle_y.at(i) - X_y.at(i));
    }

    polyfit_test(stage_angle, X, stage_angle_y, X_y);
    sin_fit_test(stage_angle, X, stage_angle_y, X_y);
    e_sin_fit_test(stage_angle, X, stage_angle_y, X_y);
    alt_sin_fit_test(stage_angle, X, stage_angle_y, X_y);
    lut_test(stage_angle, X, stage_angle_y, X_y);
    cspline_test(stage_angle, X, stage_angle_y, X_y);
    bspline_test(stage_angle, X, stage_angle_y, X_y);
}
