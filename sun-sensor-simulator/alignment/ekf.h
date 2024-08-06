#pragma once
#include <../math/matlib.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

/*
* goals:
* 1. Estimate the bias of the gyro
*/

class Filter {
public:
    Vector<double> x_hat;
    Vector<double> xHat_previous;
    Vector<double> xHat_bar;
    Vector<double> yHat;
    Matrix<double> p;
    Matrix<double> p_bar;
    Matrix<double> Q;
    Matrix<double> R;
    Matrix<double> C_m;
    Matrix<double> A;
    Matrix<double> B;
    Matrix<double> K;

private:
    Matrix<double> jacobian(Vector<double> v);
    Vector<double> predict_measured_state(Vector<double> mag_vector_model);
public:
    Filter();
    void set_quaternion_process_variance(double variance);
    void set_gyro_process_variance(double variance);
    void set_mag_vector_measurement_variance(double variance);
    void predict(double dt, Vector<double> w, Vector<double> mag_vector_model);
    void update(Vector<double> mag_vector_measured);

    Vector<double> determine_bias(std::string gyro_data_file, std::string mag_data_file, std::vector<Vector<double>> &out_mag_data, std::vector<Vector<double>> &out_gyro_data);
};
