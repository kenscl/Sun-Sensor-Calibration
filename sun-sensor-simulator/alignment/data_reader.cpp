#include "data_reader.h"
#include <ctime>
#include <fstream>
#include <sys/stat.h>
#include <vector>


Reader::Reader(std::string gyro_data_file, std::string mag_data_file, std::string sun_data_file, std::string sun_selection_file, std::string attitude_file) {
    gyro_data = read_gyro_data(gyro_data_file);
    mag_data = read_mag_data(mag_data_file);
    sun_data = read_sun_data(sun_data_file);
    sun_selection_data = read_sun_selection(sun_selection_file);
    attitude_data = read_attitude(attitude_file);
    dt = std::vector<double>();

    dt.push_back(1.);

    for (uint i = 1; i < sun_time.size(); ++i) {
        dt.push_back((sun_time.at(i) - sun_time.at(i-1)) / 1000.);
    }
}

std::vector<Vector<double>> Reader::read_gyro_data(std::string gyro_data_file) {
    std::fstream gyro_file(gyro_data_file);
    if (!gyro_file.is_open()) {
        printf("Could not open file %s\n", gyro_data_file.c_str());
        return std::vector<Vector<double>>();
    }

    std::vector<int> time;
    dt = std::vector<double>();
    gyro_data = std::vector<Vector<double>>();
    gyro_time = std::vector<int>();
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
                        gyro_time.push_back(stol(cell));
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
                return std::vector<Vector<double>>();
            } catch (const std::out_of_range& e) {
                printf("Out of range\n");
                return std::vector<Vector<double>>();
            }
        }
        gyro_data.push_back(temp);
    }
    return gyro_data;
}

std::vector<Vector<double>> Reader::read_mag_data(std::string mag_data_file){
    std::fstream mag_file(mag_data_file);
    if (!mag_file.is_open()) {
        printf("Could not open file %s\n", mag_data_file.c_str());
        return std::vector<Vector<double>>();
    }

    std::string line;
    mag_data = std::vector<Vector<double>>();
    mag_time = std::vector<int>();
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
                    case 0:
                        mag_time.push_back(stol(cell));
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
                return std::vector<Vector<double>>();
            } catch (const std::out_of_range& e) {
                printf("Out of range\n");
                return std::vector<Vector<double>>();
            }
        }
        mag_data.push_back(temp);
    }
    return mag_data;
}

std::vector<Vector<double>> Reader::read_sun_data(std::string sun_data_file) {
    std::fstream sun_file(sun_data_file);
    if (!sun_file.is_open()) { printf("Could not open file %s\n", sun_data_file.c_str());
        return std::vector<Vector<double>>();
    }

    std::string line;
    sun_data = std::vector<Vector<double>>();
    sun_time = std::vector<int>();
    // two dummy reads to start in line 3
    getline(sun_file, line);
    getline(sun_file, line);

    while (getline(sun_file, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        Vector<double> temp(3);

        int i = 0;
        while (getline(lineStream, cell, ',')) {
            try {
                switch (i) {
                    case 0:
                        sun_time.push_back(stol(cell));
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
                return std::vector<Vector<double>>();
            } catch (const std::out_of_range& e) {
                printf("Out of range\n");
                return std::vector<Vector<double>>();
            }
        }
        sun_data.push_back(temp);
    }
    return sun_data;
}

std::vector<int> Reader::read_sun_selection(std::string sun_selection_file){
    std::fstream sun_selection(sun_selection_file);
    if (!sun_selection.is_open()) {
        printf("Could not open file %s\n", sun_selection_file.c_str());
        return std::vector<int>();
    }

    sun_selection_data = std::vector<int>(); 
    sun_selection_time = std::vector<int>();
    std::string line;
    getline(sun_selection, line);
    getline(sun_selection, line);
    getline(sun_selection, line);

    while (getline(sun_selection, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        int temp;

        int i = 0;
        while (getline(lineStream, cell, ',')) {
            try {
                switch (i) {
                    case 0:
                        sun_selection_time.push_back(stol(cell));
                        break;
                    case 1:
                        temp = (stoi(cell));
                        break;
                    default:
                        break;
                }
                ++i;
            } catch (const std::invalid_argument& e) {
                printf("Invalid argument\n");
                return std::vector<int>();
            } catch (const std::out_of_range& e) {
                printf("Out of range\n");
                return std::vector<int>();
            }
        }
        sun_selection_data.push_back(temp);
    }
    return sun_selection_data;
}

std::vector<Quaternion<double>> Reader::read_attitude(std::string attitude_file) {
    std::fstream attitude(attitude_file);
    if (!attitude.is_open()) {
        printf("Could not open file %s\n", attitude_file.c_str());
        return std::vector<Quaternion<double>>();
    }
    
    attitude_data = std::vector<Quaternion<double>>();
    attitude_time = std::vector<int>();
    std::string line;
    getline(attitude, line);
    getline(attitude, line);

    while (getline(attitude, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        Quaternion<double> temp;

        int i = 0;
        while (getline(lineStream, cell, ',')) {
            try {
                switch (i) {
                    case 0:
                        attitude_time.push_back(stol(cell));
                        break;
                    case 1:
                        temp.q = (stof(cell));
                        break;
                    case 2:
                        temp.i = (stof(cell));
                        break;
                    case 3:
                        temp.j = (stof(cell));
                        break;
                    case 4:
                        temp.k = (stof(cell));
                        break;
                    default:
                        break;
                }
                ++i;
            } catch (const std::invalid_argument& e) {
                printf("Invalid argument\n");
                return std::vector<Quaternion<double>>();
            } catch (const std::out_of_range& e) {
                printf("Out of range\n");
                return std::vector<Quaternion<double>>();
            }
        }
        attitude_data.push_back(temp);
    }
    return attitude_data;
}

void Reader::reflect_data(uint count) {
    uint size = gyro_data.size() - 1;
    for (uint i = 0; i < count; ++i) {
        if (i % 2 == 0){
            for (uint j = size; j > 0; --j) {
                gyro_data.push_back(gyro_data.at(j));
                gyro_time.push_back(gyro_time.at(j));
                sun_data.push_back(sun_data.at(j));
                sun_time.push_back(sun_time.at(j));
                mag_data.push_back(mag_data.at(j));
                mag_time.push_back(mag_time.at(j));
                sun_selection_data.push_back(sun_selection_data.at(j));
                sun_selection_time.push_back(sun_selection_time.at(j));
                attitude_data.push_back(attitude_data.at(j));
                attitude_time.push_back(attitude_time.at(j));
                dt.push_back(dt.at(j));
            }
        } else {
            for (uint j = 0; j < size; ++j) {
                gyro_data.push_back(gyro_data.at(j));
                gyro_time.push_back(gyro_time.at(j));
                sun_data.push_back(sun_data.at(j));
                sun_time.push_back(sun_time.at(j));
                mag_data.push_back(mag_data.at(j));
                mag_time.push_back(mag_time.at(j));
                sun_selection_data.push_back(sun_selection_data.at(j));
                sun_selection_time.push_back(sun_selection_time.at(j));
                attitude_data.push_back(attitude_data.at(j));
                attitude_time.push_back(attitude_time.at(j));
                dt.push_back(dt.at(j));
            }
        }
    }
}
