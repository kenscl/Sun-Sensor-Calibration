#pragma once
#include "../math/matlib.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include "../lut/c-spline.h"
#include "../polynomial_fit/polynomial_fit.h"
#include "../sun_sensor_data_parser.h"
class Corrector {
private:

    std::vector<Matrix<double>> ss_SB;
    std::vector<std::string> time;
    std::vector<double> source;
    std::vector<Vector<double>> data;
    std::vector<double> angle_x_rad;
    std::vector<double> angle_y_rad;
    std::string filename;
    std::string source_fname;

public:

    Corrector(std::string filename, std::string source);

    void read_data();
    void convert_to_angle();
    void correct();
    void convert_to_vector();
    void write_data();
};





