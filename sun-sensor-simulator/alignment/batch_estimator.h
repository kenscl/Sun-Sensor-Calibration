#pragma once
#include "../math/matlib.h"
#include <vector>
class batch_estimator {
public:
    std::vector<std::vector<Vector<double>>> sun_vectors;
    std::vector<Vector<double>> mag_vectors, gyro_vectors;

    std::vector<Vector<double>> sun_refrence;
    std::vector<Vector<double>> mag_refrence, gyro_refrence;
    std::vector<int> selected_sun_vector;

    Matrix<double> dcm_SB_mag;
    Matrix<double> dcm_SB_gyro;
    std::vector<Matrix<double>> dcm_SB_ss;

    std::vector<std::vector<Vector<double>>> uncalibrated_sun_vectors;
    std::vector<Vector<double>> uncalibrated_mag_vectors;
    std::vector<Vector<double>> uncalibrated_gyro_vectors;

    std::vector<Vector<double>> cosine_error_measurements;
    std::vector<Matrix<double>> H_k;
    Vector<double> theata;

    std::vector<Matrix<double>> R_mag;
    std::vector<Matrix<double>> R_gyro;
    std::vector<std::vector<Matrix<double>>> R_sun;

   std::vector<Matrix<double>> G_mag;
   std::vector<Matrix<double>> G_gyro;
   std::vector<std::vector<Matrix<double>>> G_sun;

   std::vector<Matrix<double>> B_k;

public: 
    batch_estimator();
    void compute_uncalibrated_measurments();
    void compute_cosine_error_measurements();
    void compute_H_k();
    void comute_G_k();
    void compute_B_k();
};
