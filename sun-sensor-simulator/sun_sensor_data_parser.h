#include <cstdio>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "math/matlib.h"

inline std::vector<std::vector<double>> parse_csv(const char* file) {
    std::vector<double> stage_angle;
    std::vector<double> X;
    std::vector<double> Y;
    std::fstream data(file);
    if (!data.is_open()) {
        printf("Could not open file %s\n", file);
        return {};
    }

    std::string line;
    int line_number = 0;
    getline(data, line);
    getline(data, line);
    while (getline(data, line)) {
        std::stringstream lineStream(line);
        std::string cell;

        int i = 0;
        while (getline(lineStream, cell, ',')) {
            try {
                if (i == 1) {
                    stage_angle.push_back(stof(cell) * D2R);
                }
                if (i == 6) {
                    X.push_back(stof(cell) * D2R);
                }
                if (i == 7) {
                    Y.push_back(stof(cell) * D2R);
                }
                ++i;
            } catch (const std::invalid_argument& e) {
                return {};
            } catch (const std::out_of_range& e) {
                return {};
            }
        }

    }
    std::vector<std::vector<double>> out;
    out.push_back(stage_angle);
    out.push_back(X);
    out.push_back(Y);

    return out;
}


