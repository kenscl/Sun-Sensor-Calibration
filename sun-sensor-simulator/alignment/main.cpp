#include "sat_model.h"
#include "ekf.h"
#include <cstdio>
#include <string>
#include <vector>
#include "angle_corrector.h"
#include "data_reader.h"
#include "vasco.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
std::vector<std::vector<double>> conmpute_quat_diff(Quaternion<double> ref, std::vector<Quaternion<double>> data) {
    std::vector<Quaternion<double>> result;
    std::vector<double> q,i,j,k;
    for (uint o = 0; o < data.size(); ++o) {
        result.push_back(ref.get_diff(data.at(o)));
        i.push_back(result.at(o).to_rpy().data[0] * R2D);
        j.push_back(result.at(o).to_rpy().data[1] * R2D);
        k.push_back(result.at(o).to_rpy().data[2] * R2D);
    }
    std::vector<std::vector<double>> ret;
    ret.push_back(i);
    ret.push_back(j);
    ret.push_back(k);
    return ret;
}

std::vector<std::vector<double>> parse_vectors(std::vector<Vector<double>> v) {
    std::vector<std::vector<double>> ret(v.size());
    for (uint i = 0; i < v.at(0).size; ++i) {
        for (uint j = 0; j < v.size(); ++j) {
            ret.at(i).push_back(v.at(j).data.at(i));
        }
    }
    return ret;
}

std::vector<std::vector<double>> parse_quates(std::vector<Quaternion<double>> v) {
    std::vector<std::vector<double>> ret(4);
        for (uint j = 0; j < v.size(); ++j) {
            ret.at(0).push_back(v.at(j).q);
            ret.at(1).push_back(v.at(j).i);
            ret.at(2).push_back(v.at(j).j);
            ret.at(3).push_back(v.at(j).k);
        }
    return ret;
}

std::vector<std::vector<double>> conmpute_quat_diff(std::vector<Quaternion<double>> ref, std::vector<Quaternion<double>> data) {
    std::vector<Quaternion<double>> result;
    std::vector<double> q,i,j,k;
    for (uint o = 0; o < data.size(); ++o) {
        result.push_back(ref.data()[o].get_diff(data.at(o)));
        i.push_back(result.at(o).to_rpy().data[0] * R2D);
        j.push_back(result.at(o).to_rpy().data[1] * R2D);
        k.push_back(result.at(o).to_rpy().data[2] * R2D);
    }
    std::vector<std::vector<double>> ret;
    ret.push_back(i);
    ret.push_back(j);
    ret.push_back(k);
    return ret;
}

void plot_vasco(Vasco &vasco, std::string name, Quaternion<double> mag_misaignemnt, std::vector<Quaternion<double>> ss_misalignments, std::vector<Quaternion<double>> attitude = {}) {
    
    std::vector<std::vector<double>> gt = parse_quates(attitude);
    std::vector<std::vector<double>> estimates = parse_quates(vasco.attitude_estimates);
    plt::cla();
    plt::clf();
    plt::title("ATTIUDE ESTIMATION ERROR");
    plt::xlabel("TIME [s]");
    plt::ylabel("ERROR VALUE [deg]");
    plt::named_plot("ground truth q", vasco.elapsed_time, gt[0]);
    plt::named_plot("ground truth i", vasco.elapsed_time, gt[1]);
    plt::named_plot("ground truth j", vasco.elapsed_time, gt[2]);
    plt::named_plot("ground truth k", vasco.elapsed_time, gt[3]);
    plt::named_plot("estimate q", vasco.elapsed_time, estimates[0]);
    plt::named_plot("estimate i", vasco.elapsed_time, estimates[1]);
    plt::named_plot("estimate j", vasco.elapsed_time, estimates[2]);
    plt::named_plot("estimate k", vasco.elapsed_time, estimates[3]);
    plt::legend();
    plt::save("../plots/" + name + "att_err.eps");

    std::vector<std::vector<double>> mag_error = conmpute_quat_diff(Quaternion<double>(1,0,0,0), vasco.mag_misalignment_estimates);
    plt::cla();
    plt::clf();
    plt::title(name + " MAGNETOMETER");
    plt::xlabel("TIME [s]");
    plt::ylabel("ERROR VALUE [deg]");
    plt::named_plot("x", vasco.elapsed_time, mag_error[0]);
    plt::named_plot("y", vasco.elapsed_time, mag_error[1]);
    plt::named_plot("z", vasco.elapsed_time, mag_error[2]);
    plt::save("../plots/" + name + "_mag_misalignment.eps");

    plt::cla();
    plt::clf();
    std::vector<std::vector<double>> ss_XP_error = conmpute_quat_diff(Quaternion<double>(1,0,0,0), vasco.sun_sensor_xp_misalignment_estimates);
    plt::title(name + " SUN SENSOR XP MISALIGNMENT");
    plt::xlabel("TIME [s]");
    plt::ylabel("ERROR VALUE [deg]");
    plt::named_plot("x", vasco.elapsed_time, ss_XP_error[0]);
    plt::named_plot("y", vasco.elapsed_time, ss_XP_error[1]);
    plt::named_plot("z", vasco.elapsed_time, ss_XP_error[2]);
    plt::save("../plots/" + name + "_ss_xp_misalignment.eps");


    plt::cla();
    plt::clf();
    std::vector<std::vector<double>> ss_XM_error = conmpute_quat_diff(Quaternion<double>(1,0,0,0), vasco.sun_sensor_xm_misalignment_estimates);
    plt::title(name + " SUN SENSOR XM");
    plt::xlabel("TIME [s]");
    plt::ylabel("ERROR VALUE [deg]");
    plt::named_plot("x", vasco.elapsed_time, ss_XM_error[0]);
    plt::named_plot("y", vasco.elapsed_time, ss_XM_error[1]);
    plt::named_plot("z", vasco.elapsed_time, ss_XM_error[2]);
    plt::save("../plots/" + name + "_ss_xm_misalignment.eps");

    plt::cla();
    plt::clf();
    std::vector<std::vector<double>> ss_YP_error = conmpute_quat_diff(Quaternion<double>(1,0,0,0), vasco.sun_sensor_yp_misalignment_estimates);
    plt::title(name + " SUN SENSOR YP");
    plt::xlabel("TIME [s]");
    plt::ylabel("ERROR VALUE [deg]");
    plt::named_plot("x", vasco.elapsed_time, ss_YP_error[0]);
    plt::named_plot("y", vasco.elapsed_time, ss_YP_error[1]);
    plt::named_plot("z", vasco.elapsed_time, ss_YP_error[2]);
    plt::save("../plots/" + name + "_ss_yp_misalignment.eps");

    plt::cla();
    plt::clf();
    std::vector<std::vector<double>> ss_YM_error = conmpute_quat_diff(Quaternion<double>(1,0,0,0), vasco.sun_sensor_ym_misalignment_estimates);
    plt::title(name + " SUN SENSOR YM");
    plt::xlabel("TIME [s]");
    plt::ylabel("ERROR VALUE [deg]");
    plt::named_plot("x", vasco.elapsed_time, ss_YM_error[0]);
    plt::named_plot("y", vasco.elapsed_time, ss_YM_error[1]);
    plt::named_plot("z", vasco.elapsed_time, ss_YM_error[2]);
    plt::save("../plots/" + name + "_ss_ym_misalignment.eps");

    plt::cla();
    plt::clf();
    std::vector<std::vector<double>> ss_ZP_error = conmpute_quat_diff(Quaternion<double>(1,0,0,0), vasco.sun_sensor_zp_misalignment_estimates);
    plt::title(name + " SUN SENSOR ZP");
    plt::xlabel("TIME [s]");
    plt::ylabel("ERROR VALUE [deg]");
    plt::named_plot("x", vasco.elapsed_time, ss_ZP_error[0]);
    plt::named_plot("y", vasco.elapsed_time, ss_ZP_error[1]);
    plt::named_plot("z", vasco.elapsed_time, ss_ZP_error[2]);
    plt::save("../plots/" + name + "_ss_zp_misalignment.eps");

    plt::cla();
    plt::clf();
    std::vector<std::vector<double>> ss_ZM_error = conmpute_quat_diff(Quaternion<double>(1,0,0,0), vasco.sun_sensor_zm_misalignment_estimates);
    plt::title(name + " SUN SENSOR ZM");
    plt::xlabel("TIME [s]");
    plt::ylabel("ERROR VALUE [deg]");
    plt::named_plot("x", vasco.elapsed_time, ss_ZM_error[0]);
    plt::named_plot("y", vasco.elapsed_time, ss_ZM_error[1]);
    plt::named_plot("z", vasco.elapsed_time, ss_ZM_error[2]);
    plt::save("../plots/" + name + "_ss_zm_misalignment.eps");

    
    std::vector<std::vector<double>> cov_entry(48);
    for (uint i = 0; i < vasco.cov_matrix_diag.size(); ++i) {
        for (uint j = 0; j < 48; ++j) {
            cov_entry.at(j).push_back(vasco.cov_matrix_diag.at(i).data.at(j));
        }
    }

    plt::cla();
    plt::clf();
    plt::title("COVARIANCE MATRIX");
    plt::xlabel("TIME [s]");
    plt::ylabel("COVARIANCE VALUE");
    for (uint i = 0; i < 48; ++i) {
        plt::named_plot("entry " + std::to_string(i), vasco.elapsed_time, cov_entry.at(i));
    }
    plt::save("../plots/" + name + "_covariance_matrix.eps");
    plt::show();

    std::vector<std::vector<double>> residuals = parse_vectors(vasco.residuals);

    plt::cla();
    plt::clf();
    plt::title(name + " Residuals");
    plt::xlabel("TIME [s]");
    plt::ylabel("RESIDUAL VALUE []");
    for (uint i = 0; i < 21; ++i) {
        plt::named_plot("entry " + std::to_string(i), vasco.elapsed_time, residuals.at(i));
    }
    plt::save("../plots/" + name + "residuals.eps");
}

void determine_misalignment_z_x() {
    Reader data("../../full_sat_data/z_x/Gyro-data-as-joinbyfield-2024-07-26 18_27_32.csv",
                "../../full_sat_data/z_x/Magvector-data-as-joinbyfield-2024-07-26 18_27_48.csv", 
                "../../full_sat_data/z_x/Sun Vector Body Selected-data-as-joinbyfield-2024-07-26 18_27_59.csv_corrected.csv",
                "../../full_sat_data/z_x/Source-data-2024-07-29 11_54_52.csv",
                "../../full_sat_data/z_x/Quaternion-data-as-joinbyfield-2024-07-26 18_27_20.csv");
    printf("sizes: %lu, %lu, %lu, %lu, %lu\n", data.gyro_data.size(), data.mag_data.size(), data.sun_data.size(), data.sun_selection_data.size(), data.attitude_data.size());
    data.reflect_data(200);

    std::vector<Vector<double>> sun_refrence, mag_refrence;
    Vector<double> sun_model(3), mag_model(3);
    mag_model = data.mag_data[0];
    sun_model = data.sun_data[data.sun_data.size()/2];

    for (uint i = 0; i < data.gyro_data.size(); ++i) {
        sun_refrence.push_back(sun_model);
        mag_refrence.push_back(mag_model);
        data.attitude_data[i] = data.attitude_data[i].normalize();
    }
    data.mag_data[0].print();
    data.mag_data[data.mag_data.size() - 1].print();

    std::vector<Quaternion<double>> q_ss;
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));

    Vasco vasco;
    vasco.run(data.dt, data.gyro_data, mag_refrence, data.mag_data, sun_refrence, data.sun_data, data.sun_selection_data, 2, data.attitude_data);
    vasco.full_feedback(-1);
    vasco.print(-1);
    
    plot_vasco(vasco, "MISALIGNMENT Z-X -", Quaternion<double>(1,0,0,0), q_ss);
}

void determine_misalignment_z_y() {
    Reader data("../../full_sat_data/z_y/Gyro-data-as-joinbyfield-2024-07-26 18_50_51.csv",
                "../../full_sat_data/z_y/Magvector-data-as-joinbyfield-2024-07-26 18_51_01.csv", 
                "../../full_sat_data/z_y/Sun Vector Body Selected-data-as-joinbyfield-2024-07-26 18_51_09.csv_corrected.csv",
                "../../full_sat_data/z_y/Source-data-2024-07-29 11_52_31.csv",
                "../../full_sat_data/z_y/Quaternion-data-as-joinbyfield-2024-07-26 18_50_31.csv");
    printf("sizes: %lu, %lu, %lu, %lu, %lu\n", data.gyro_data.size(), data.mag_data.size(), data.sun_data.size(), data.sun_selection_data.size(), data.attitude_data.size());
    data.reflect_data(200);
    std::vector<Vector<double>> sun_refrence, mag_refrence;
    Vector<double> sun_model(3), mag_model(3);
    mag_model = data.mag_data[0];
    sun_model = data.sun_data[data.sun_data.size()/2];

    for (uint i = 0; i < data.gyro_data.size(); ++i) {
        sun_refrence.push_back(sun_model);
        mag_refrence.push_back(mag_model);
        data.attitude_data[i] = data.attitude_data[i].normalize();
    }
    data.mag_data[0].print();
    data.mag_data[data.mag_data.size() - 1].print();

    std::vector<Quaternion<double>> q_ss;
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));

    Vasco vasco;
    vasco.run(data.dt, data.gyro_data, mag_refrence, data.mag_data, sun_refrence, data.sun_data, data.sun_selection_data, 2, data.attitude_data);
    vasco.full_feedback(-1);
    vasco.print(-1);
    
    plot_vasco(vasco, "MISALIGNMENT Z-Y -", Quaternion<double>(1,0,0,0), q_ss);

}

void determine_misalignment_xzx() {
    Reader data("../../full_sat_data/xzx/Gyro-data-as-joinbyfield-2024-07-26 17_40_24.csv",
                "../../full_sat_data/xzx/Magvector-data-as-joinbyfield-2024-07-26 17_40_31.csv", 
                "../../full_sat_data/xzx/Sun Vector Body Selected-data-as-joinbyfield-2024-07-26 17_40_42.csv",
                "../../full_sat_data/xzx/Source-data-2024-07-29 11_55_31.csv",
                "../../full_sat_data/xzx/Quaternion-data-as-joinbyfield-2024-07-26 17_40_13.csv");
    printf("sizes: %lu, %lu, %lu, %lu, %lu\n", data.gyro_data.size(), data.mag_data.size(), data.sun_data.size(), data.sun_selection_data.size(), data.attitude_data.size());
    data.reflect_data(200);

    std::vector<Vector<double>> sun_refrence, mag_refrence;
    Vector<double> sun_model(3), mag_model(3);
    mag_model = data.mag_data[0];
    sun_model = data.sun_data[data.sun_data.size()/2];

    for (uint i = 0; i < data.gyro_data.size(); ++i) {
        sun_refrence.push_back(sun_model);
        mag_refrence.push_back(mag_model);
        data.attitude_data[i] = data.attitude_data[i].normalize();
    }
    data.mag_data[0].print();
    data.mag_data[data.mag_data.size() - 1].print();

    std::vector<Quaternion<double>> q_ss;
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));

    Vasco vasco;
    vasco.run(data.dt, data.gyro_data, mag_refrence, data.mag_data, sun_refrence, data.sun_data, data.sun_selection_data, 2, data.attitude_data);
    vasco.full_feedback(-1);
    vasco.print(-1);
    
    plot_vasco(vasco, "MISALIGNMENT Y-Z-Y -", Quaternion<double>(1,0,0,0), q_ss);
}

void determine_misalignment_yzy() {
    Reader data("../../full_sat_data/yzy/Gyro-data-as-joinbyfield-2024-07-26 17_07_11.csv",
                "../../full_sat_data/yzy/Magvector-data-as-joinbyfield-2024-07-26 17_07_20.csv", 
                "../../full_sat_data/yzy/Sun Vector Body Selected-data-as-joinbyfield-2024-07-26 17_07_30.csv",
                "../../full_sat_data/yzy/Source-data-2024-07-29 11_57_01.csv",
                "../../full_sat_data/yzy/Quaternion-data-as-joinbyfield-2024-07-26 17_07_02.csv");
    printf("sizes: %lu, %lu, %lu, %lu, %lu\n", data.gyro_data.size(), data.mag_data.size(), data.sun_data.size(), data.sun_selection_data.size(), data.attitude_data.size());
    data.reflect_data(200);

    std::vector<Vector<double>> sun_refrence, mag_refrence;
    Vector<double> sun_model(3), mag_model(3);
    mag_model = data.mag_data[0];
    sun_model = data.sun_data[data.sun_data.size()/2];

    for (uint i = 0; i < data.gyro_data.size(); ++i) {
        sun_refrence.push_back(sun_model);
        mag_refrence.push_back(mag_model);
        data.attitude_data[i] = data.attitude_data[i].normalize();
    }
    data.mag_data[0].print();
    data.mag_data[data.mag_data.size() - 1].print();

    std::vector<Quaternion<double>> q_ss;
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));

    Vasco vasco;
    vasco.run(data.dt, data.gyro_data, mag_refrence, data.mag_data, sun_refrence, data.sun_data, data.sun_selection_data, 2, data.attitude_data);
    vasco.full_feedback(-1);
    vasco.print(-1);
    
    plot_vasco(vasco, "MISALIGNMENT X-Z-X -", Quaternion<double>(1,0,0,0), q_ss);
}

void determine_misalignment_xyxy() {
    Reader data("../../full_sat_data/yxyx/Gyro-data-as-joinbyfield-2024-07-26 18_40_43.csv",
                "../../full_sat_data/yxyx/Magvector-data-as-joinbyfield-2024-07-26 18_40_51.csv",
                "../../full_sat_data/yxyx/Sun Vector Body Selected-data-as-joinbyfield-2024-07-26 18_41_00.csv",
                "../../full_sat_data/yxyx/Sun_selection_generated.csv",
                "../../full_sat_data/yxyx/Quaternion-data-as-joinbyfield-2024-07-26 18_40_34.csv");
    printf("sizes: %lu, %lu, %lu, %lu, %lu\n", data.gyro_data.size(), data.mag_data.size(), data.sun_data.size(), data.sun_selection_data.size(), data.attitude_data.size());
    data.reflect_data(200);

    std::vector<Vector<double>> sun_refrence, mag_refrence;
    Vector<double> sun_model(3), mag_model(3);
    mag_model = data.mag_data[0];
    sun_model = data.sun_data[data.sun_data.size()/2];

    for (uint i = 0; i < data.gyro_data.size(); ++i) {
        sun_refrence.push_back(sun_model);
        mag_refrence.push_back(mag_model);
        data.attitude_data[i] = data.attitude_data[i].normalize();
    }
    data.mag_data[0].print();
    data.mag_data[data.mag_data.size() - 1].print();

    std::vector<Quaternion<double>> q_ss;
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));

    Vasco vasco;
    vasco.run(data.dt, data.gyro_data, mag_refrence, data.mag_data, sun_refrence, data.sun_data, data.sun_selection_data, 2, data.attitude_data);
    vasco.full_feedback(-1);
    vasco.print(-1);
    
    plot_vasco(vasco, "MISALIGNMENT Y-X-Y-X -", Quaternion<double>(1,0,0,0), q_ss);
}

void sim_no_misalignment() {

    Matrix<double> mag_misalignment, gyro_misalignments;
    std::vector<Matrix<double>> ss_misaligments;
    std::vector<Matrix<double>> ss_alignments;
    std::vector<Quaternion<double>> q_ss;
    ss_alignments.push_back(Matrix<double>().Identity(3));
    ss_alignments.push_back(Matrix<double>().Identity(3));
    ss_alignments.push_back(Matrix<double>().Identity(3));
    ss_alignments.push_back(Matrix<double>().Identity(3));
    ss_alignments.push_back(Matrix<double>().Identity(3));
    ss_alignments.push_back(Matrix<double>().Identity(3));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));

    Vector<double> sun_vector(3);
    sun_vector.data.at(0) = 5;
    sun_vector.data.at(1) = 1;
    sun_vector = sun_vector.normalize();

    Vector<double> mag_vector(3);
    mag_vector.data.at(0) = 3;
    mag_vector.data.at(1) = 0;
    mag_vector.data.at(2) = 1;
    mag_vector = mag_vector.normalize();

    Vector<double> axis(3);
    axis.data.at(1) = 1;

    Quaternion<double> q_mag( 0* D2R, axis);
    q_mag = q_mag.normalize();
    mag_misalignment = q_mag.to_rotation_matrix();
    mag_misalignment.print();

    Sat_Model sat_model(mag_misalignment, Matrix<double>().Identity(3), ss_alignments, sun_vector, mag_vector);
    std::vector<Vector<double>> sun_meas;
    std::vector<Vector<double>> sun_vectors, sun_model;
    std::vector<Vector<double>> mag_meas, mag_model;
    std::vector<Vector<double>> gyro_meas;
    std::vector<int> selected_sun_vector;
    std::vector<double> dt;
    std::vector<Quaternion<double>> attitude;


    Vector<double> rotation(3);
    for (int p = 0; p < 20; ++p) {
        for (int i = 0; i < 3; ++i) {
            switch (i) {
                case 0:
                    rotation.data.at(0) = 1;
                    rotation.data.at(1) = 0;
                    rotation.data.at(2) = 0;
                    break; 
                case 1:
                    rotation.data.at(0) = 0;
                    rotation.data.at(1) = 1;
                    rotation.data.at(2) = 0;
                    break;
                case 2:
                    rotation.data.at(0) = 0;
                    rotation.data.at(1) = 0;
                    rotation.data.at(2) = 1;
                    break;
                case 3:
                    rotation.data.at(0) = 1;
                    rotation.data.at(1) = 0;
                    rotation.data.at(2) = 0;
                    break;
                case 4:
                    rotation.data.at(0) = 0;
                    rotation.data.at(1) = 1;
                    rotation.data.at(2) = 0;
                    break;
                case 5:
                    rotation.data.at(0) = 0;
                    rotation.data.at(1) = 0;
                    rotation.data.at(2) = 1;
                    break;
            }

            Quaternion<double> gyro_integrated(1,0,0,0);

            for (double j = 0; j < 360 * D2R; j+=5*D2R) {
                for (int u = 0; u < 5; ++u) {
                    bool has_value = false;
                    sun_meas = sat_model.get_sun_vectors();
                    for (uint k = 0; k < sun_meas.size(); ++k) {
                        if (sun_meas.at(k).norm() > 0.01) {
                            sun_vectors.push_back(sun_meas.at(k));
                            sun_model.push_back(sat_model.sun_vector);
                            selected_sun_vector.push_back(k + 1);
                            has_value = true;
                        }
                    }

                    if (!has_value) {
                        sun_vectors.push_back(Vector<double>(3));
                        sun_model.push_back(Vector<double>(3)); 
                        selected_sun_vector.push_back(8);
                    }
                    gyro_meas.push_back(sat_model.get_gyro_vector());
                    Quaternion<double> gyro_q(gyro_meas.at(gyro_meas.size() - 1));
                    mag_meas.push_back(sat_model.get_mag_vector());
                    mag_model.push_back(sat_model.mag_vector);
                    dt.push_back(1);
                    attitude.push_back(sat_model.attitude);

                }
                sat_model.rotate(5 * D2R, rotation);

            }

        }

    }

    Vasco vasco;
    vasco.run(dt, gyro_meas, mag_model, mag_meas, sun_model, sun_vectors, selected_sun_vector, 15, attitude);
    plot_vasco(vasco, "SIMULATION NO MISALIGNMENT -", q_mag, q_ss);
    vasco.print(-1);

}

void sim_mag_misalignment() {

    Matrix<double> mag_misalignment, gyro_misalignments;
    std::vector<Matrix<double>> ss_misaligments;
    std::vector<Matrix<double>> ss_alignments;
    std::vector<Quaternion<double>> q_ss;
    ss_alignments.push_back(Matrix<double>().Identity(3));
    ss_alignments.push_back(Matrix<double>().Identity(3));
    ss_alignments.push_back(Matrix<double>().Identity(3));
    ss_alignments.push_back(Matrix<double>().Identity(3));
    ss_alignments.push_back(Matrix<double>().Identity(3));
    ss_alignments.push_back(Matrix<double>().Identity(3));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));
    q_ss.push_back(Quaternion<double>(Matrix<double>().Identity(3)));

    Vector<double> sun_vector(3);
    sun_vector.data.at(0) = 5;
    sun_vector.data.at(1) = 1;
    sun_vector = sun_vector.normalize();

    Vector<double> mag_vector(3);
    mag_vector.data.at(0) = 3;
    mag_vector.data.at(1) = 0;
    mag_vector.data.at(2) = 1;
    mag_vector = mag_vector.normalize();

    Vector<double> axis(3);
    axis.data.at(1) = 1;

    Quaternion<double> q_mag( 3* D2R, axis);
    q_mag = q_mag.normalize();
    mag_misalignment = q_mag.to_rotation_matrix();
    mag_misalignment.print();

    Sat_Model sat_model(mag_misalignment, Matrix<double>().Identity(3), ss_alignments, sun_vector, mag_vector);
    std::vector<Vector<double>> sun_meas;
    std::vector<Vector<double>> sun_vectors, sun_model;
    std::vector<Vector<double>> mag_meas, mag_model;
    std::vector<Vector<double>> gyro_meas;
    std::vector<int> selected_sun_vector;
    std::vector<double> dt;
    std::vector<Quaternion<double>> attitude;


    Vector<double> rotation(3);
    for (int p = 0; p < 20; ++p) {
        for (int i = 0; i < 3; ++i) {
            switch (i) {
                case 0:
                    rotation.data.at(0) = 1;
                    rotation.data.at(1) = 0;
                    rotation.data.at(2) = 0;
                    break; 
                case 1:
                    rotation.data.at(0) = 0;
                    rotation.data.at(1) = 1;
                    rotation.data.at(2) = 0;
                    break;
                case 2:
                    rotation.data.at(0) = 0;
                    rotation.data.at(1) = 0;
                    rotation.data.at(2) = 1;
                    break;
                case 3:
                    rotation.data.at(0) = 1;
                    rotation.data.at(1) = 0;
                    rotation.data.at(2) = 0;
                    break;
                case 4:
                    rotation.data.at(0) = 0;
                    rotation.data.at(1) = 1;
                    rotation.data.at(2) = 0;
                    break;
                case 5:
                    rotation.data.at(0) = 0;
                    rotation.data.at(1) = 0;
                    rotation.data.at(2) = 1;
                    break;
            }

            Quaternion<double> gyro_integrated(1,0,0,0);

            for (double j = 0; j < 360 * D2R; j+=5*D2R) {
                for (int u = 0; u < 5; ++u) {
                    bool has_value = false;
                    sun_meas = sat_model.get_sun_vectors();
                    for (uint k = 0; k < sun_meas.size(); ++k) {
                        if (sun_meas.at(k).norm() > 0.01) {
                            sun_vectors.push_back(sun_meas.at(k));
                            sun_model.push_back(sat_model.sun_vector);
                            selected_sun_vector.push_back(k + 1);
                            has_value = true;
                        }
                    }

                    if (!has_value) {
                        sun_vectors.push_back(Vector<double>(3));
                        sun_model.push_back(Vector<double>(3)); 
                        selected_sun_vector.push_back(8);
                    }
                    gyro_meas.push_back(sat_model.get_gyro_vector());
                    Quaternion<double> gyro_q(gyro_meas.at(gyro_meas.size() - 1));
                    mag_meas.push_back(sat_model.get_mag_vector());
                    mag_model.push_back(sat_model.mag_vector);
                    dt.push_back(1);
                    attitude.push_back(sat_model.attitude);

                }
                sat_model.rotate(5 * D2R, rotation);

            }

        }

    }

    Vasco vasco;
    vasco.run(dt, gyro_meas, mag_model, mag_meas, sun_model, sun_vectors, selected_sun_vector, 1, attitude);

    plot_vasco(vasco, "SIMULATION MAGNETOMETER MISALIGNMENT -", q_mag, q_ss);
    vasco.print(-1);
}

void sim_all_misalignment() {

    Matrix<double> mag_misalignment, gyro_misalignments;
    std::vector<Matrix<double>> ss_misaligments;
    std::vector<Matrix<double>> ss_alignments;
    std::vector<Quaternion<double>> q_ss;

    Vector<double> sun_vector(3);
    sun_vector.data.at(0) = 5;
    sun_vector.data.at(1) = 1;
    sun_vector = sun_vector.normalize();

    Vector<double> mag_vector(3);
    mag_vector.data.at(0) = 3;
    mag_vector.data.at(1) = 0;
    mag_vector.data.at(2) = 1;
    mag_vector = mag_vector.normalize();

    Vector<double> axis(3), xp_axis(3), xm_axis(3), yp_axis(3), ym_axis(3), zp_axis(3), zm_axis(3);
    axis.data.at(1) = 1;
    xp_axis.data.at(0) = 1;

    xm_axis.data.at(0) = -1;
    xm_axis.data.at(1) = -1;

    yp_axis.data.at(1) = 1;
    yp_axis.data.at(2) = 1;

    ym_axis.data.at(1) = -1;
    ym_axis.data.at(2) = 1;

    zp_axis.data.at(2) = 1;
    zp_axis.data.at(0) = 1;

    zm_axis.data.at(2) = -1;
    zm_axis.data.at(0) = 1;

    Quaternion<double> q_mag( 3* D2R, axis);
    Quaternion<double> q_ss_xp( 3* D2R, xp_axis);
    Quaternion<double> q_ss_xm( 3* D2R, xm_axis);
    Quaternion<double> q_ss_yp( 3* D2R, yp_axis);
    Quaternion<double> q_ss_ym( 3* D2R, ym_axis);
    Quaternion<double> q_ss_zp( 3* D2R, zp_axis);
    Quaternion<double> q_ss_zm( 3* D2R, zm_axis);

    q_mag = q_mag.normalize();
    q_ss_xp = q_ss_xp.normalize();
    q_ss_xm = q_ss_xm.normalize();
    q_ss_yp = q_ss_yp.normalize();
    q_ss_ym = q_ss_ym.normalize();
    q_ss_zp = q_ss_zp.normalize();
    q_ss_zm = q_ss_zm.normalize();

    mag_misalignment = q_mag.to_rotation_matrix();
    ss_misaligments.push_back(q_ss_xp.to_rotation_matrix());
    ss_misaligments.push_back(q_ss_xm.to_rotation_matrix());
    ss_misaligments.push_back(q_ss_yp.to_rotation_matrix());
    ss_misaligments.push_back(q_ss_ym.to_rotation_matrix());
    ss_misaligments.push_back(q_ss_zp.to_rotation_matrix());
    ss_misaligments.push_back(q_ss_zm.to_rotation_matrix());

    q_ss.push_back(q_ss_xp);
    q_ss.push_back(q_ss_xm);
    q_ss.push_back(q_ss_yp);
    q_ss.push_back(q_ss_ym);
    q_ss.push_back(q_ss_zp);
    q_ss.push_back(q_ss_zm);

    mag_misalignment.print();

    Sat_Model sat_model(mag_misalignment, Matrix<double>().Identity(3), ss_misaligments, sun_vector, mag_vector);
    std::vector<Vector<double>> sun_meas;
    std::vector<Vector<double>> sun_vectors, sun_model;
    std::vector<Vector<double>> mag_meas, mag_model;
    std::vector<Vector<double>> gyro_meas;
    std::vector<int> selected_sun_vector;
    std::vector<double> dt;
    std::vector<Quaternion<double>> attitude;


    Vector<double> rotation(3);
    for (int p = 0; p < 10; ++p) {
        for (int i = 0; i < 3; ++i) {
            switch (i) {
                case 0:
                    rotation.data.at(0) = 1;
                    rotation.data.at(1) = 0;
                    rotation.data.at(2) = 0;
                    break; 
                case 1:
                    rotation.data.at(0) = 0;
                    rotation.data.at(1) = 1;
                    rotation.data.at(2) = 0;
                    break;
                case 2:
                    rotation.data.at(0) = 0;
                    rotation.data.at(1) = 0;
                    rotation.data.at(2) = 1;
                    break;
                case 3:
                    rotation.data.at(0) = 1;
                    rotation.data.at(1) = 0;
                    rotation.data.at(2) = 0;
                    break;
                case 4:
                    rotation.data.at(0) = 0;
                    rotation.data.at(1) = 1;
                    rotation.data.at(2) = 0;
                    break;
                case 5:
                    rotation.data.at(0) = 0;
                    rotation.data.at(1) = 0;
                    rotation.data.at(2) = 1;
                    break;
            }

            Quaternion<double> gyro_integrated(1,0,0,0);

            for (double j = 0; j < 360 * D2R; j+=5*D2R) {
                for (int u = 0; u < 5; ++u) {
                    bool has_value = false;
                    sun_meas = sat_model.get_sun_vectors();
                    for (uint k = 0; k < sun_meas.size(); ++k) {
                        if (sun_meas.at(k).norm() > 0.01) {
                            sun_vectors.push_back(sun_meas.at(k));
                            sun_model.push_back(sat_model.sun_vector);
                            selected_sun_vector.push_back(k + 1);
                            has_value = true;
                        }
                    }

                    if (!has_value) {
                        sun_vectors.push_back(Vector<double>(3));
                        sun_model.push_back(Vector<double>(3)); 
                        selected_sun_vector.push_back(8);
                    }
                    gyro_meas.push_back(sat_model.get_gyro_vector());
                    Quaternion<double> gyro_q(gyro_meas.at(gyro_meas.size() - 1));
                    mag_meas.push_back(sat_model.get_mag_vector());
                    mag_model.push_back(sat_model.mag_vector);
                    dt.push_back(1);
                    attitude.push_back(sat_model.attitude);

                }
                sat_model.rotate(5 * D2R, rotation);

            }

        }

    }

    Vasco vasco;
    vasco.run(dt, gyro_meas, mag_model, mag_meas, sun_model, sun_vectors, selected_sun_vector, 1, attitude);

    plot_vasco(vasco, "ALL SENSOR MISALIGNMENT -", q_mag, q_ss);
    vasco.print(-1);
}

int main () {
    sim_all_misalignment();
}

