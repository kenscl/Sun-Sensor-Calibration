#include "batch_estimator.h"

batch_estimator::batch_estimator() {
    dcm_SB_mag = Matrix<double>(3,3).Identity(3);
    dcm_SB_gyro = Matrix<double>(3,3).Identity(3);
    dcm_SB_ss = std::vector<Matrix<double>>();

    Matrix<double> ss_xp_dcm_SB(3,3);
    ss_xp_dcm_SB.data.at(1).at(0) = 1;
    ss_xp_dcm_SB.data.at(2).at(1) = 1;
    ss_xp_dcm_SB.data.at(0).at(2) = 1;
    dcm_SB_ss.push_back(ss_xp_dcm_SB);
    
    Matrix<double> ss_xm_dcm_SB(3,3);
    ss_xm_dcm_SB.data.at(1).at(0) = -1;
    ss_xm_dcm_SB.data.at(2).at(1) = 1;
    ss_xm_dcm_SB.data.at(0).at(2) = -1;
    dcm_SB_ss.push_back(ss_xm_dcm_SB);

    Matrix<double> ss_yp_dcm_SB(3,3);
    ss_yp_dcm_SB.data.at(0).at(0) = -1;
    ss_yp_dcm_SB.data.at(2).at(1) = 1;
    ss_yp_dcm_SB.data.at(1).at(2) = 1;
    dcm_SB_ss.push_back(ss_yp_dcm_SB);

    Matrix<double> ss_ym_dcm_SB(3,3);
    ss_ym_dcm_SB.data.at(0).at(0) = 1;
    ss_ym_dcm_SB.data.at(2).at(1) = 1;
    ss_ym_dcm_SB.data.at(1).at(2) = -1;
    dcm_SB_ss.push_back(ss_ym_dcm_SB);

    Matrix<double> ss_zp_dcm_SB(3,3);
    ss_zp_dcm_SB.data.at(0).at(0) = -1;
    ss_zp_dcm_SB.data.at(1).at(1) = -1;
    ss_zp_dcm_SB.data.at(2).at(2) = 1;
    dcm_SB_ss.push_back(ss_zp_dcm_SB);

    Matrix<double> ss_zm_dcm_SB(3,3);
    ss_zm_dcm_SB.data.at(0).at(0) = 1;
    ss_zm_dcm_SB.data.at(1).at(1) = -1;
    ss_zm_dcm_SB.data.at(2).at(2) = -1;
    dcm_SB_ss.push_back(ss_zm_dcm_SB);
}

void batch_estimator::compute_uncalibrated_measurments() {
    uncalibrated_sun_vectors = std::vector<std::vector<Vector<double>>>();
    uncalibrated_mag_vectors = std::vector<Vector<double>>();
    uncalibrated_gyro_vectors = std::vector<Vector<double>>();

    for (uint i = 0; i < sun_vectors.size(); i++) {
        std::vector<Vector<double>> sun_vector;
        for (uint j = 0; j < dcm_SB_ss.size(); j++) {
            sun_vector.push_back(dcm_SB_ss.at(j) * sun_vectors.at(i).at(j));
        }
        uncalibrated_sun_vectors.push_back(sun_vector);
    }

    for (uint i = 0; i < mag_vectors.size(); i++) {
        uncalibrated_mag_vectors.push_back(dcm_SB_mag * mag_vectors.at(i));
    }

    for (uint i = 0; i < gyro_vectors.size(); i++) {
        uncalibrated_gyro_vectors.push_back(dcm_SB_gyro * gyro_vectors.at(i));
    }
}

void batch_estimator::compute_cosine_error_measurements() {
    cosine_error_measurements = std::vector<Vector<double>>();
    for (uint i = 0; i < sun_refrence.size(); i++) {
        Vector<double> cosine_error(3);
        printf("if you got here you forgot to integrate the gyro ! \n");
        cosine_error.data[0] = uncalibrated_sun_vectors[i][selected_sun_vector[i]] * uncalibrated_mag_vectors[i] - sun_refrence[i] * mag_refrence[i];
        cosine_error.data[1] = uncalibrated_sun_vectors[i][selected_sun_vector[i]] * uncalibrated_gyro_vectors[i] - sun_refrence[i] * gyro_vectors[i];
        cosine_error.data[2] = uncalibrated_mag_vectors[i] * uncalibrated_gyro_vectors[i] - mag_refrence[i] * gyro_vectors[i];
        cosine_error_measurements.push_back(cosine_error);
    }
}

void batch_estimator::compute_H_k(){
    Matrix<double> H(3,9);
    for (uint i = 0; i < cosine_error_measurements.size(); i++) {
        Vector<double> a = cosine_error_measurements.at(i);
        H = a.to_matrix() * theata.transpose() * (1 / (theata.transpose() * theata).data[0]);
    }
    H_k.push_back(H);
}

void batch_estimator::comute_G_k(){
    for (uint i = 0; i < cosine_error_measurements.size(); i++) {
        G_gyro.push_back(R_gyro[i].cholesky());
        G_mag.push_back(R_mag[i].cholesky());
        for (uint j = 0; j < R_sun[i].size(); j++) {
            G_sun[i].push_back(R_sun[i][j].cholesky());
        }
    }
}

void batch_estimator::compute_B_k() {
    for (uint i = 0; i < cosine_error_measurements.size(); i++) {
        Matrix<double> b_k_temp = Matrix<double>(3, 9);
        Matrix<double> b_k = Matrix<double>(3, 9);
        Matrix<double> b_gyro = gyro_vectors[i].transpose() * G_gyro[i];
        Matrix<double> b_mag = mag_vectors[i].transpose() * G_mag[i];
        Matrix<double> b_sun = sun_vectors[i][selected_sun_vector[i]].transpose() * G_sun[i][selected_sun_vector[i]];
        concatenate_axis_0(b_gyro, b_mag, b_k_temp);
        concatenate_axis_0(b_k_temp, b_sun, b_k);
        B_k.push_back(b_k);
    }
}
