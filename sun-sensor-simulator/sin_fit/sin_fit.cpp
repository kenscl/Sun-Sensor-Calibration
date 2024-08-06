#include "sin_fit.h"
#include <cstdio>

Matrix<double> Sin_Fit::Jacobian(int64_t degree, const std::vector<double>& gt_data) {
    Matrix<double> jacobian(gt_data.size(), 2 * degree); 

    for (size_t j = 0; j < gt_data.size(); j++) {
        double x = gt_data.at(j);
        for (int i = 0; i < degree; i++) {
            double amp_param = this->parameters.data.at(2 * i);  
            double freq_param = this->parameters.data.at(2 * i + 1);

            jacobian.data.at(j).at(2 * i) = sin(freq_param * x); // dR/d(b1)
            jacobian.data.at(j).at(2 * i + 1) = amp_param * x * cos(freq_param * x); // dR/d(b2)
        }
    }
    return jacobian;
}

Sin_Fit::Sin_Fit(int64_t degree, int64_t max_steps, double bound, std::vector<double> gt_data, std::vector<double> measurement) {
    this->degree = degree;

    double ak = .1;
    double bk = .001;
    Vector<double> params = Vector<double>((int) 2 * degree);
    this->parameters = params;
    for (int i = 0; i < params.size; i++) {
        if (i %2 == 0) this->parameters.data.at(i) = 1.;
        else this->parameters.data.at(i) = 0.0000000000000000000001;
    } 

    /*
     * 1. Inital guess
     * 2. calculate the residuals
     * 3. S = sum r^2
     */

    double S = MAXFLOAT;
    int steps = 0;
    while (S > bound) {
        // end regression if to many steps are done
        if (steps > max_steps) break;
        // calculating the residuals
        Vector<double> residuals(measurement.size());
        for (uint j = 0; j < measurement.size(); j++) {
            residuals.data.at(j) = measurement.at(j) - this->sf_at(gt_data.at(j));
        }

        // calculating S
        S = 0;
        for (uint j = 0; j < residuals.size; j++) {
            S += residuals.data.at(j) * residuals.data.at(j);
        }

        Matrix<double> jacobian = this->Jacobian(degree, gt_data);

        Matrix<double> I = Matrix<double>::Identity(jacobian.cols);

        Matrix<double> JtJ = jacobian.transpose() * jacobian;
        Matrix<double> JtJ_reg = JtJ + I * bk;

        Matrix<double> inv = JtJ_reg.inverse();

        Vector<double> delta_params = inv * jacobian.transpose() * residuals * ak;
        this->parameters = this->parameters + delta_params;
        double S_new = 0;

        for (uint j = 0; j < measurement.size(); j++) {
            double new_residual = measurement.at(j) - this->sf_at(gt_data.at(j));
            S_new += new_residual * new_residual;
        }

        if (S_new < S) {
            bk /= 10; 
        } else {
            bk *= 10; 
        }

        // end regression if no notable improvement happens
        double delta_s = abs(S_new - S);
        if (delta_s < bound / 1000) break;

        ++steps; 
    }
    printf("residual at end: %f after %d steps \n", S, steps);

}

double Sin_Fit::sf_at(double x) {
    double sum = 0;
    for (int k = 0; k < this->degree; ++k) {
        double b0 = this->parameters.data.at(2 * k);
        double b1 = this->parameters.data.at(2 * k + 1);
        sum += b0 * sin(b1 * x);
    }
    return sum;
}


double Sin_Fit::sf_at_derivative(double x) {
    double sum = 0;
    for (int k = 0; k < this->degree; ++k) {
        double b0 = this->parameters.data.at(2 * k);
        double b1 = this->parameters.data.at(2 * k + 1);
        sum += b0 * b1 * cos(b1 * x);
    }
    return sum;
}
std::vector<double> Sin_Fit::sf_calc(std::vector<double> x) {
    std::vector<double> ret;
    for (uint i = 0; i < x.size(); i++) {
        ret.push_back(this->sf_at(x.at(i)));
    }
    return ret;
}
double Sin_Fit::at(double x, double bound) {
    double delta = MAXFLOAT;
    double last_x = x;

    for (int i = 0; i < 20; i++ ){
    //while (delta > bound) {
        double f = this->sf_at(last_x);
        double df = this->sf_at_derivative(last_x);
        double new_x;
        new_x = last_x - (f - x) / df;
        delta = abs(new_x - last_x);
        last_x = new_x;
    }
    
    return last_x;
}

std::vector<double> Sin_Fit::calc(std::vector<double> x, double bound) {
    std::vector<double> ret;
    for (uint i = 0; i < x.size(); i++) {
        ret.push_back(this->at(x.at(i), bound));
    }
    return ret;
}
