#pragma once
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

#include "../polynomial_fit/polynomial_fit.h"
#include "../sin_fit/sin_fit.h"
#include "../sin_fit/extended_sin_fit.h"
#include "../sin_fit/extended_sin_fit_alt.h"
#include "../lut/lut.h"
#include "../lut/c-spline.h"
#include "../lut/b-spline.h"



#include <chrono>
namespace timer = std::chrono;
#include <vector>

void split_vector_odd_even(std::vector<double> init, std::vector<double>& even, std::vector<double>& odd);

void polyfit_test(std::vector<double> gt_x, std::vector<double> measurement_x, std::vector<double> gt_y, std::vector<double> measurement_y);
void sin_fit_test(std::vector<double> gt_x, std::vector<double> measurement_x, std::vector<double> gt_y, std::vector<double> measurement_y);
void e_sin_fit_test(std::vector<double> gt_x, std::vector<double> measurement_x, std::vector<double> gt_y, std::vector<double> measurement_y);
void alt_sin_fit_test(std::vector<double> gt_x, std::vector<double> measurement_x, std::vector<double> gt_y, std::vector<double> measurement_y);
void lut_test(std::vector<double> gt_x, std::vector<double> measurement_x, std::vector<double> gt_y, std::vector<double> measurement_y);
void cspline_test(std::vector<double> gt_x, std::vector<double> measurement_x, std::vector<double> gt_y, std::vector<double> measurement_y);
void bspline_test(std::vector<double> gt_x, std::vector<double> measurement_x, std::vector<double> gt_y, std::vector<double> measurement_y);

