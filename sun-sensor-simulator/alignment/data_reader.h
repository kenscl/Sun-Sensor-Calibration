#include <string>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include "../math/matlib.h"

class Reader {
public:
    std::vector<Vector<double>> gyro_data, mag_data, sun_data;
    std::vector<Quaternion<double>> attitude_data;
    std::vector<double> dt;
    std::vector<int> sun_selection_data;
    std::vector<int> gyro_time, mag_time, sun_time, sun_selection_time, attitude_time; 

    Reader(std::string gyro_data_file, std::string mag_data_file, std::string sun_data_file, std::string sun_selection_file, std::string attitude_file);
    std::vector<Vector<double>> read_gyro_data(std::string gyro_data_file);
    std::vector<Vector<double>> read_mag_data(std::string mag_data_file);
    std::vector<Vector<double>> read_sun_data(std::string sun_data_file);
    std::vector<int> read_sun_selection(std::string sun_selection_file);
    std::vector<Quaternion<double>> read_attitude(std::string attitude_file);
    
    void reflect_data(uint count);
};
