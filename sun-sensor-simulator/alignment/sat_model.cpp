#include "sat_model.h"
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <system_error>
Sat_Model::Sat_Model(Matrix<double> mag_misalignment, Matrix<double> gyro_misalignments, std::vector<Matrix<double>> ss_misaligments, Vector<double> sun_vector, Vector<double> mag_vector) {
    this->attitude = Quaternion<double>(1, 0, 0, 0);
    this->gyro_measurment = Vector<double>(3);
    this->mag_misalignment = mag_misalignment;
    this->gyro_misalignments = gyro_misalignments;
    this->ss_misaligments = ss_misaligments;

    this->sun_vector = sun_vector;
    this->mag_vector = mag_vector;

    // init Sensor to Body dcms
    
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

std::vector<Vector<double>> Sat_Model::get_sun_vectors() {
    std::vector<Vector<double>> ret;
    bool found = false;
    for (uint i = 0; i < this->dcm_SB_ss.size(); i++) {
        Vector<double> normal, normal_S(3), measured;
        normal_S.data.at(2) = 1;

        Quaternion<double> rotated_reference(0, sun_vector.data[0], sun_vector.data[1], sun_vector.data[2]);
        Vector<double> sHat_N(3);

        // Rotate the reference vector by the quaternion
        rotated_reference = attitude.conjugate() * rotated_reference * attitude;
        sHat_N.data[0] = rotated_reference.i;
        sHat_N.data[1] = rotated_reference.j;
        sHat_N.data[2] = rotated_reference.k;

        Vector<double> sHat_S = (dcm_SB_ss.at(i)).inverse()  * sHat_N;

        measured = ss_misaligments.at(i) * dcm_SB_ss.at(i) * sHat_S;

        double angle = sHat_S.calculate_angle(normal_S);
        
        if (angle > 45 * D2R) {
            measured.data.at(0) = 0.;
            measured.data.at(1) = 0.;
            measured.data.at(2) = 0.;
            ret.push_back(measured);
        } else if (angle < -45 * D2R) {
            measured.data.at(0) = 0.;
            measured.data.at(1) = 0.;
            measured.data.at(2) = 0.;
            ret.push_back(measured);
        } else if (!found) {
            ret.push_back(measured.normalize());
            found = true;
        } else {
            measured.data.at(0) = 0.;
            measured.data.at(1) = 0.;
            measured.data.at(2) = 0.;
            ret.push_back(measured);
        }
    }
    return ret;
}

Vector<double> Sat_Model::get_mag_vector() {
    Quaternion<double> rotated_reference(0, mag_vector.data[0], mag_vector.data[1], mag_vector.data[2]);
    Vector<double> mHat_N(3);

    // Rotate the reference vector by the quaternion
    rotated_reference = attitude.conjugate() * rotated_reference * attitude;

    // Extract the vector part of the rotated quaternion
    mHat_N.data[0] = rotated_reference.i;
    mHat_N.data[1] = rotated_reference.j;
    mHat_N.data[2] = rotated_reference.k;
    // Transform the vector using the inverse of the DCM and then back
    Vector<double> mHat_S = (dcm_SB_mag).inverse() * mHat_N;
    Vector<double> result = mag_misalignment * dcm_SB_mag * mHat_S;

    return result.normalize();
}

Vector<double> Sat_Model::get_gyro_vector() {
    Vector<double> result = gyro_measurment;
    gyro_measurment = Vector<double>(3);
    return result;
}

void Sat_Model::rotate(double angle, Vector<double> axis) {
    // change attitude for sensors
    Quaternion<double> rotation(angle, axis);
    Quaternion<double> last_attitude = attitude;
    attitude = rotation * attitude;
    attitude = attitude.normalize();

    // determihne gyro measurement
    gyro_measurment = axis * angle;
}
