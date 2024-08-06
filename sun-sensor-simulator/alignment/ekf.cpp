#include "ekf.h"
#include <cstdio>
#include <string>
#include <vector>

Filter::Filter() {
    x_hat = Vector<double>(7);
    x_hat.data[0] = 0.5;
    x_hat.data[1] = 0.5;
    x_hat.data[2] = 0.5;
    x_hat.data[3] = 0.5;
    x_hat.data[4] = 0.;
    x_hat.data[5] = 0.;
    x_hat.data[6] = 0.;
    yHat = Vector<double>(3);
    p = Matrix<double>(7,7);
    p = p.I() * 0.1;
    Q = Matrix<double>(7,7);
    Q = Q.I();
    R = Matrix<double>(3,3);
    R = R.I();
    double quaternion_process_variance = 0.001;

    double gyro_bias_variance = 1;
    set_quaternion_process_variance(quaternion_process_variance);
    set_gyro_process_variance(gyro_bias_variance);
    double mag_vector_variance = 0.1;
    set_mag_vector_measurement_variance(mag_vector_variance);
    A = Matrix<double>(7,7);
    B = Matrix<double>(7,6);
    C_m = Matrix<double>(3,7);
}

void Filter::set_quaternion_process_variance(double variance) {
    Vector<double> diagonal = Q.diag();
    diagonal.data[0] = variance;
    diagonal.data[1] = variance;
    diagonal.data[2] = variance;
    diagonal.data[3] = variance;
    diagonal.data[4] = diagonal.data[4];
    diagonal.data[5] = diagonal.data[5];
    diagonal.data[6] = diagonal.data[6];
    Q.diag(diagonal);
}

void Filter::set_gyro_process_variance(double variance) {
    Vector<double> diagonal = Q.diag();
    diagonal.data[0] = diagonal.data[0];
    diagonal.data[1] = diagonal.data[1];
    diagonal.data[2] = diagonal.data[2];
    diagonal.data[3] = diagonal.data[3];
    diagonal.data[4] = variance;
    diagonal.data[5] = variance;
    diagonal.data[6] = variance;
    Q.diag(diagonal);
}

void Filter::set_mag_vector_measurement_variance(double variance) {
    Vector<double> diagonal = R.diag();
    diagonal.data[0] = variance;
    diagonal.data[1] = variance;
    diagonal.data[2] = variance;
    R.diag(diagonal);
}

Matrix<double> Filter::jacobian(Vector<double> v) {
    double x = v.data[0];
    double y = v.data[1];
    double z = v.data[2];
    double q0 = xHat_previous.data[0];
    double q1 = xHat_previous.data[1];
    double q2 = xHat_previous.data[2];
    double q3 = xHat_previous.data[3];

    Matrix<double> jacobian(3,4);
    jacobian.data[0][0] = x * q0 - z * q2 + y * q3;
    jacobian.data[0][1] = x * q1 + y * q2 + z * q3;
    jacobian.data[0][2] = -z * q0 + y * q1 - x * q2;
    jacobian.data[0][3] =  y * q0 + z * q1 - x * q3;

    jacobian.data[1][0] = y * q0 + z * q1 - x * q3;
    jacobian.data[1][1] = z * q0 - y * q1 + x * q2;
    jacobian.data[1][2] = x * q1 + y * q2 + z * q3;
    jacobian.data[1][3] = -x * q0 + z * q2 - y * q3;

    jacobian.data[2][0] = z * q0 - y * q1 + x * q2;
    jacobian.data[2][1] = -y * q0 - z * q1 + x * q3;
    jacobian.data[2][2] = x * q0 - z * q2 + y * q3;
    jacobian.data[2][3] = x * q1 + y * q2 + z * q3;

    jacobian = jacobian * 2;
    return jacobian;

}
Vector<double> Filter::predict_measured_state(Vector<double> mag_vector_model) {
    double q0 = xHat_bar.data[0];
    double q1 = xHat_bar.data[1];
    double q2 = xHat_bar.data[2];
    double q3 = xHat_bar.data[3];
    Matrix<double> rot_mat(3,3);

    rot_mat.data[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    rot_mat.data[0][1] = 2 * (q1 * q2 - q0 * q3);
    rot_mat.data[0][2] = 2 * (q1 * q3 + q0 * q2);
    rot_mat.data[1][0] = 2 * (q1 * q2 + q0 * q3);
    rot_mat.data[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    rot_mat.data[1][2] = 2 * (q2 * q3 - q0 * q1);
    rot_mat.data[2][0] = 2 * (q1 * q3 - q0 * q2);
    rot_mat.data[2][1] = 2 * (q2 * q3 + q0 * q1);
    rot_mat.data[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
    rot_mat = rot_mat.transpose();

    Matrix<double> hPrime_mag = this->jacobian(mag_vector_model);
    Vector<double> mag_Bar = rot_mat * mag_vector_model;
    Matrix<double> zeros(3,3);
    concatenate_axis_1<double>(hPrime_mag, zeros, C_m);
    //C_m = hPrime_mag;
    return mag_Bar;
}

void Filter::predict(double dt, Vector<double> w, Vector<double> mag_vector_model) {
    Quaternion<double> q(x_hat.data.at(0), x_hat.data.at(1), x_hat.data.at(2), x_hat.data.at(3));
    Matrix<double> s_q(4,3);
    s_q.data.at(0).at(0) = -q.i;
    s_q.data.at(0).at(1) = -q.j;
    s_q.data.at(0).at(2) = -q.k;

    s_q.data.at(1).at(0) = q.q;
    s_q.data.at(1).at(1) = -q.k;
    s_q.data.at(1).at(2) = q.j;
    
    s_q.data.at(2).at(0) = q.k;
    s_q.data.at(2).at(1) = q.q;
    s_q.data.at(2).at(2) = -q.i;

    s_q.data.at(3).at(0) = -q.j;
    s_q.data.at(3).at(1) = q.i;
    s_q.data.at(3).at(2) = q.q;

    Matrix<double> I4(4,4);
    I4 = I4.I();
    Matrix<double> tmp1(4,7);
    Matrix<double> Sq_time_dependent_neg = Matrix<double>(s_q * (-dt /2));
    concatenate_axis_1<double>(I4,Sq_time_dependent_neg, tmp1);

    Matrix<double> I3(3,3);
    I3 = I3.I();
    Matrix<double> zeroes_3_4(3,4);
    Matrix<double> tmp2(3,7);
    concatenate_axis_1<double>(zeroes_3_4, I3, tmp2);

    concatenate_axis_0<double>(tmp1,tmp2, A);

    Matrix<double> zeroes_3_3(3,3);
    Matrix<double> Sq_time_dependent = Matrix<double>(s_q * (dt / 2));
    concatenate_axis_0<double>(Sq_time_dependent, zeroes_3_3, B);
    Vector<double> measurement(6);
    measurement.data[0] = w.data[0];
    measurement.data[1] = w.data[1];
    measurement.data[2] = w.data[2];
    xHat_bar = A * x_hat + B * measurement;

    Quaternion<double> q_temp(xHat_bar.data[0], xHat_bar.data[1], xHat_bar.data[2], xHat_bar.data[3]);
    q_temp = q_temp.normalize();
    xHat_bar.data[0] = q_temp.q;
    xHat_bar.data[1] = q_temp.i;
    xHat_bar.data[2] = q_temp.j;
    xHat_bar.data[3] = q_temp.k;
    xHat_previous = xHat_bar;

    yHat = this->predict_measured_state(mag_vector_model);
    p_bar = A * p * A.transpose() + Q;
}

void Filter::update(Vector<double> mag_vector_measured) {
    Matrix<double> tmp1 = C_m * p_bar * C_m.transpose() + R;
    tmp1 = tmp1.inverse();
    K = p_bar * C_m.transpose() * tmp1;


    x_hat = xHat_bar + K * (mag_vector_measured - yHat);
    Quaternion<double> q(x_hat.data.at(0), x_hat.data.at(1), x_hat.data.at(2), x_hat.data.at(3));
    q = q.normalize();
    x_hat.data[0] = q.q;
    x_hat.data[1] = q.i;
    x_hat.data[2] = q.j;
    x_hat.data[3] = q.k;

    Matrix<double> I7(7,7);
    I7 = I7.I();
    p = (I7 - K * C_m) * p_bar;
}

Vector<double> Filter::determine_bias(std::string gyro_data_file, std::string mag_data_file, std::vector<Vector<double>> &out_mag_data, std::vector<Vector<double>> &out_gyro_data) {
    std::fstream gyro_file(gyro_data_file);
    if (!gyro_file.is_open()) {
        printf("Could not open file %s\n", gyro_data_file.c_str());
        return Vector<double>(3);
    }

    std::fstream mag_file(mag_data_file);
    if (!mag_file.is_open()) {
        printf("Could not open file %s\n", mag_data_file.c_str());
        return Vector<double>(3);
    }


    std::vector<int> time;
    std::vector<Vector<double>> gyro_data;
    std::string line;

    // two dummy reads to start in line 3
    getline(gyro_file, line);
    getline(gyro_file, line);
    while (getline(gyro_file, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        Vector<double> temp(3);

        int i = 0;
        while (getline(lineStream, cell, ',')) {
            try {
                switch (i) {
                    case 0:
                        time.push_back(stol(cell));
                        break;
                    case 1:
                        temp.data[0] = (stof(cell));
                        break;
                    case 2:
                        temp.data[1] = (stof(cell));
                        break;
                    case 3:
                        temp.data[2] = (stof(cell));
                        break;
                    default:
                        break;
                }
                ++i;
            } catch (const std::invalid_argument& e) {
                printf("Invalid argument\n");
                return Vector<double>(3);
            } catch (const std::out_of_range& e) {
                printf("Out of range\n");
                return Vector<double>(3);
            }
        }
        gyro_data.push_back(temp);
    }

    std::vector<Vector<double>> mag_data;
    // two dummy reads to start in line 3
    getline(mag_file, line);
    getline(mag_file, line);

    while (getline(mag_file, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        Vector<double> temp(3);

        int i = 0;
        while (getline(lineStream, cell, ',')) {
            try {
                switch (i) {
                    case 4:
                        temp.data[0] = (stof(cell));
                        break;
                    case 5:
                        temp.data[1] = (stof(cell));
                        break;
                    case 6:
                        temp.data[2] = (stof(cell));
                        break;
                    default:
                        break;
                }
                ++i;
            } catch (const std::invalid_argument& e) {
                printf("Invalid argument\n");
                return Vector<double>(3);
            } catch (const std::out_of_range& e) {
                printf("Out of range\n");
                return Vector<double>(3);
            }
        }
        mag_data.push_back(temp);
    }

    std::vector<double> dt;

    dt.push_back(1.);

    for (uint i = 1; i < time.size(); ++i) {
        dt.push_back((time.at(i) - time.at(i-1)) / 1000.);
    }

    //determine bias

    /*
    * refrence magnetic field is that of Würzburg in NEO, but actual magnetic field dosnt matter
    */
    Vector<double> mag_vector_wuerzburg(3);
    mag_vector_wuerzburg.data[0] = 20.6675;
    mag_vector_wuerzburg.data[1] = -0.861;
    mag_vector_wuerzburg.data[2] = 43.6435;
    mag_vector_wuerzburg = mag_vector_wuerzburg.normalize();

    for (uint i = 0; i < mag_data.size(); ++i) {
        Vector<double> w = gyro_data.at(i);
        Vector<double> mag = mag_data.at(i);
        this->predict(dt.at(i), w, mag_vector_wuerzburg);
        this->update(mag.normalize());
    }

    Vector<double> bias(3);
    bias.data[0] = x_hat.data[4];
    bias.data[1] = x_hat.data[5];
    bias.data[2] = x_hat.data[6];

    for (uint i = 0; i < mag_data.size(); ++i) {
        mag_data.at(i).data[0] -= bias.data[0];
    }

    out_mag_data = mag_data;
    out_gyro_data = gyro_data;
    return bias;

}
