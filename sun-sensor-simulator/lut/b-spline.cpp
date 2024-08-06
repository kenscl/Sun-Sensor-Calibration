#include "b-spline.h"
#include <strings.h>
#include <vector>


BSpline::BSpline(std::vector<double> gt_data, std::vector<double> measurement, int num_points) {
    std::vector<double> x, y;
    this->basis_matrix = Matrix<double>(4,4);
    basis_matrix.data[0][0] = 1./6.;
    basis_matrix.data[0][1] = 4./6.;
    basis_matrix.data[0][2] = 1./6.;

    basis_matrix.data[1][0] = -3./6.;
    basis_matrix.data[1][2] = 3./6.;

    basis_matrix.data[2][0] = 3./6.;
    basis_matrix.data[2][1] = -6./6.;
    basis_matrix.data[2][2] = 3./6.;


    basis_matrix.data[3][0] = -1./6.;
    basis_matrix.data[3][1] = 3./6.;
    basis_matrix.data[3][2] = -3./6.;
    basis_matrix.data[3][3] = 1./6.;

    for (int i = 0; i < gt_data.size(); ++i) {
        y.push_back(gt_data.at(i) - measurement.at(i));
        x.push_back(measurement.at(i));
    }
    this->x = select_evenly_spaced_points(x, num_points);
    this->y = select_evenly_spaced_points(y, num_points);
}

std::vector<double> BSpline::get_value(double x) {
    int index = this->y.size() - 5;
    for (int i = 0; i < this->y.size() - 4; ++i) {
        if (x <= this->x.at(0)) {
            index = 0;
            break;
        }
        if (this->x.at(i) <= x && x <= this->x.at(i+1)) {
            index = i;
            break;
        }
    }


    double t = (x- this->x.at(index)) / (this->x.at(index + 4) - this->x.at(index));
    
    Vector<double> x_vec = Vector<double>(4);
    for (int i = 0; i < 4; ++i) {
        x_vec.data[i] = pow(t , i);
    }
    Matrix<double> x_mat = x_vec.transpose();

    Vector<double> x_values(4), y_values(4);
    for (int i = 0; i < 4; ++i) {
        x_values.data[i] = this->x.at(index + i);
        y_values.data[i] = this->y.at(index + i);
    }

    Vector<double> result_x = x_mat * this->basis_matrix * x_values;
    Vector<double> result_y = x_mat * this->basis_matrix * y_values;

    printf("val %f, t %f, x: %f, y: %f\n",x,t,  result_x.data[0], result_y.data[0]);
    return {result_x.data[0], result_y.data[0]};
}


std::vector<std::vector<double>> BSpline::get_value(std::vector<double> x) {
    std::vector<double> x_values, y_values;
    for (int i = 0; i < x.size(); ++i) {
        std::vector<double> result = this->get_value(x.at(i));
        x_values.push_back(result.at(0));
        y_values.push_back(result.at(1));

    }
    return {x_values, y_values};
}


double BSpline::at(double x) {
    return x + this->get_value(x).at(1);
}

std::vector<double> BSpline::calc(std::vector<double> x) {
    std::vector<double> ret;
    for (uint i = 0; i < x.size(); i++) {
        ret.push_back(this->at(x.at(i)));
    }
    return ret;
}
