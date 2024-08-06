#include "c-spline.h"
#include <cstdio>
#include <vector>

CSpline::CSpline(std::vector<double> gt_data, std::vector<double> measurement, int num_points) {
    std::vector<double> x, y, mid_x, mid_y;


    y.push_back(gt_data.at(0) - measurement.at(0));
    x.push_back(measurement.at(0));
    for (int i = 0; i < gt_data.size(); ++i) {
        if (gt_data.at(i) < - 47 * D2R || gt_data.at(i) > 47 * D2R){
            continue;
        }
        else {
            mid_y.push_back(gt_data.at(i) - measurement.at(i));
            mid_x.push_back(measurement.at(i));
        }
    }
    mid_x = select_evenly_spaced_points(mid_x, num_points - 2);
    mid_y = select_evenly_spaced_points(mid_y, num_points - 2);
    x.insert(x.end(), mid_x.begin(), mid_x.end());
    y.insert(y.end(), mid_y.begin(), mid_y.end());
    y.push_back(gt_data.at(gt_data.size() - 1) - measurement.at(gt_data.size() - 1));
    x.push_back(measurement.at(gt_data.size() - 1));


    this->coefficients = generate_coefficients(x, y);
}


std::vector<std::vector<double>> CSpline::generate_coefficients(std::vector<double> x, std::vector<double> y) {
    // caclulate h
    std::vector<double> h;
    for (int i = 0; i < x.size() - 1; ++i) {
        h.push_back(x.at(i+1) - x.at(i));
    }
    // generate matrix A
    Matrix<double> A(x.size(), x.size());
    A.data[0][0] = 1;
    A.data[x.size() - 1][x.size() - 1] = 1;

    for (int i = 1; i < x.size() - 1; ++i) {
        A.data[i][i-1] = h.at(i-1);
        A.data[i][i] = 2 * (h.at(i-1) + h.at(i));
        A.data[i][i+1] = h.at(i);
    }

    // generate vector
    Vector<double> b(x.size());
    b.data[0] = 0;
    b.data[x.size() - 1] = 0;
    for (uint i = 1; i < b.size - 1; ++i) {
        b.data[i] = 3 / h.at(i) * (y.at(i+1) - y.at(i)) - 3 / h.at(i-1) * (y.at(i) - y.at(i-1));
    }

    Vector<double> c = A.inverse() * b;
    // generate coefficients
    Vector<double> a(x.size()), b_vec(x.size()), d(x.size());
    for (int i = 0; i < x.size(); ++i) {
        a.data[i] = y [i];
    }

    for ( int i = 0; i < x.size() - 1; ++i) {
        b_vec.data[i] = 1 / h.at(i) * (a.data[i+1] - a.data[i]) - h.at(i) / 3 * (2 * c.data[i] + c.data[i+1]); 
        d.data[i] = (c.data[i+1] - c.data[i]) / (3 * h.at(i));
    }

    std::vector<std::vector<double>> coefficients;
    for (uint i = 0; i < x.size(); ++i) {
        std::vector<double> temp;
        temp.push_back(x.at(i));
        temp.push_back(a.data[i]);
        temp.push_back(b_vec.data[i]);
        temp.push_back(c.data[i]);
        temp.push_back(d.data[i]);
        coefficients.push_back(temp);
    }
    return coefficients;
}


double CSpline::get_value(double x) {
    for (uint i = 0; i < this->coefficients.size() - 1 ; ++i) {
        if (x >= this->coefficients.at(i).at(0) && x <= this->coefficients.at(i+1).at(0)) {
            double x_i = this->coefficients.at(i).at(0);
            double a = this->coefficients.at(i).at(1);
            double b = this->coefficients.at(i).at(2);
            double c = this->coefficients.at(i).at(3);
            double d = this->coefficients.at(i).at(4);
            return a + b * (x - x_i) + c * (x - x_i) * (x - x_i) + d * (x - x_i) * (x - x_i) * (x - x_i);
        }
    }
    double x_i = this->coefficients.at(this->coefficients.size() - 1).at(0);
    double a = this->coefficients.at(this->coefficients.size() - 1).at(1);
    double b = this->coefficients.at(this->coefficients.size() - 1).at(2);
    double c = this->coefficients.at(this->coefficients.size() - 1).at(3);
    double d = this->coefficients.at(this->coefficients.size() - 1).at(4);
    return a + b * (x - x_i) + c * (x - x_i) * (x - x_i) + d * (x - x_i) * (x - x_i) * (x - x_i);

    return 0;
}

std::vector<double> CSpline::get_value(std::vector<double> x) {
std::vector<double> ret;
    for (uint i = 0; i < x.size(); ++i) {
        ret.push_back(this->get_value(x.at(i)));
    }
    return ret;
}

double CSpline::at(double x) {
    return x + this->get_value(x);
}

std::vector<double> CSpline::calc(std::vector<double> x) {
    std::vector<double> ret;
    for (uint i = 0; i < x.size(); i++) {
        ret.push_back(this->at(x.at(i)));
    }
    return ret;
}
