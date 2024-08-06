#include "angle_corrector.h"
#include <cmath>
#include <cstdio>

Corrector::Corrector(std::string filename, std::string source) {
    this->filename = filename;
    this->source_fname = source;

    // define matricies 

    ss_SB = std::vector<Matrix<double>>();

    Matrix<double> ss_xp_dcm_SB(3,3);
    ss_xp_dcm_SB.data.at(1).at(0) = 1;
    ss_xp_dcm_SB.data.at(2).at(1) = 1;
    ss_xp_dcm_SB.data.at(0).at(2) = 1;
    ss_SB.push_back(ss_xp_dcm_SB);
    
    Matrix<double> ss_xm_dcm_SB(3,3);
    ss_xm_dcm_SB.data.at(1).at(0) = -1;
    ss_xm_dcm_SB.data.at(2).at(1) = 1;
    ss_xm_dcm_SB.data.at(0).at(2) = -1;
    ss_SB.push_back(ss_xm_dcm_SB);

    Matrix<double> ss_yp_dcm_SB(3,3);
    ss_yp_dcm_SB.data.at(0).at(0) = -1;
    ss_yp_dcm_SB.data.at(2).at(1) = 1;
    ss_yp_dcm_SB.data.at(1).at(2) = 1;
    ss_SB.push_back(ss_yp_dcm_SB);

    Matrix<double> ss_ym_dcm_SB(3,3);
    ss_ym_dcm_SB.data.at(0).at(0) = 1;
    ss_ym_dcm_SB.data.at(2).at(1) = 1;
    ss_ym_dcm_SB.data.at(1).at(2) = -1;
    ss_SB.push_back(ss_ym_dcm_SB);

    Matrix<double> ss_zp_dcm_SB(3,3);
    ss_zp_dcm_SB.data.at(0).at(0) = -1;
    ss_zp_dcm_SB.data.at(1).at(1) = -1;
    ss_zp_dcm_SB.data.at(2).at(2) = 1;
    ss_SB.push_back(ss_zp_dcm_SB);

    Matrix<double> ss_zm_dcm_SB(3,3);
    ss_zm_dcm_SB.data.at(0).at(0) = 1;
    ss_zm_dcm_SB.data.at(1).at(1) = -1;
    ss_zm_dcm_SB.data.at(2).at(2) = -1;
    ss_SB.push_back(ss_zm_dcm_SB);

    read_data();
    convert_to_angle();
    correct();
    convert_to_vector();
    write_data();
}

void Corrector::read_data() {
    // assuming combined data
    std::fstream file(filename);
    std::fstream source(source_fname);
    if (!file.is_open()) {
        printf("Could not open file %s\n", filename.c_str());
        return;
    }
    if (!source.is_open()) {
        printf("Could not open file %s\n", source_fname.c_str());
        return;
    }

    std::string line, line_source;

    getline(file, line);
    while (getline(file, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        Vector<double> temp(3);

        int i = 0;
        while (getline(lineStream, cell, ',')) {
            try {
                switch (i) {
                    case 0:
                        time.push_back(cell);
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
                return;
            } catch (const std::out_of_range& e) {
                return;
            }
        }
        data.push_back(temp);
    }

    getline(source, line_source);
    getline(source, line_source);
    while (getline(source, line_source)) {
        std::stringstream lineStream(line_source);
        std::string cell;


        int i = 0;
        while (getline(lineStream, cell, ',')) {
            try {
                switch (i) {
                    case 0:
                        break;
                    case 1:
                        this->source.push_back(stof(cell) - 1);
                        break;
                    default:
                        break;
                }
                ++i;
            } catch (const std::invalid_argument& e) {
                return;
            } catch (const std::out_of_range& e) {
                return;
            }
        }
    }

}

void Corrector::convert_to_angle() {
    //first transform all data to local frame
    for (int i = 0; i < data.size(); ++i) {
        if (source[i] < 6) { 
            Matrix<double> ss_SB_inv = ss_SB[source[i]].inverse();
            data[i] = ss_SB_inv * data[i].normalize();
        }
    }

    // then we calculate the angles
    for (int i = 0; i < data.size(); i++) {
        angle_x_rad.push_back(M_PI / 2 - atan2(data[i].data[2], data[i].data[0]));
        angle_y_rad.push_back(M_PI / 2 - atan2(data[i].data[2], data[i].data[1]));
    }

}

void Corrector::correct() {
    std::vector<std::vector<double>> data, data_valid;
    data = parse_csv("../../Sunsensor_Data/28.08.2021 17.20.54 ADPD2140_X.txt");
    std::vector<double> stage_angle = data.at(0);
    std::vector<double> X = data.at(1);

    data_valid = parse_csv("../../Sunsensor_Data/28.08.2021 16.37.25 ADPD2140_Y.txt");
    std::vector<double> stage_angle_y = data_valid.at(0);
    std::vector<double> X_y= data_valid.at(2);

    CSpline spline_x(stage_angle, X, 76);
    CSpline spline_y(stage_angle_y, X_y, 48);
    
    angle_x_rad = spline_x.calc(angle_x_rad);
    angle_y_rad = spline_y.calc(angle_y_rad);
}

void Corrector::convert_to_vector() {
    // first convert to vector
    for (int i = 0; i < angle_x_rad.size(); ++i) {
        double angle_x = M_PI / 2 - this->angle_x_rad[i];
        double angle_y = M_PI / 2 - this->angle_y_rad[i];

        double cos_y_pow = cos(angle_y)*cos(angle_y);
        double sin_x_pow = sin(angle_x)*sin(angle_x);
        double sin_y_pow = sin(angle_y)*sin(angle_y);
        double temp_sqrt = sqrt(cos_y_pow*sin_x_pow+sin_y_pow);
        
        if (data[i].norm() > 0.00000000001) {
            data[i].data[0] = cos(angle_x)*sin(angle_y)/temp_sqrt;
            data[i].data[1] = cos(angle_y)*sin(angle_x)/temp_sqrt;
            data[i].data[2] = sin(angle_x)*sin(angle_y)/temp_sqrt;
        }

        if (temp_sqrt != 0) {
            //printf("X: %f, Y: %f, Z: %f\n", data[i].data[0], data[i].data[1], data[i].data[2]);
            data[i] = data[i].normalize();
        } else {
            //printf("angle_x: %f, angle_y: %f\n", angle_x_rad[i], angle_y_rad[i]);
            data[i].data[0] = 0;
            data[i].data[1] = 0;
            data[i].data[2] = 0;
        }
    }

    // then convert to global frame
    for (int i = 0; i < data.size(); ++i) {
        if (source[i] < 6) { 
            data[i] = ss_SB[source[i]] * data[i];
        }
        if (source[i] >= 6) {
            data[i] = Vector<double>(3);
        }
    }
} 

void Corrector::write_data() {
    std::ofstream file;
    file.open(filename + "_corrected.csv");
    file << "Time, X, Y, Z\n";
    for (int i = 0; i < data.size(); i++) {
        file << time[i] << "," << data[i].data[0] << "," << data[i].data[1] << "," << data[i].data[2] << "\n";
    }
    file.close();
}
