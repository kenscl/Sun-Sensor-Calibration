#include "extended_sin_fit.h"

Matrix<double> Extended_sin_fit::Jacobian(const std::vector<double> &gt_data) {
    double b1 = this->parameters.data.at(1);
    double b2 = this->parameters.data.at(2);
    double b3 = this->parameters.data.at(3);
    double b4 = this->parameters.data.at(4);

    Matrix<double> jacobian(gt_data.size(), 5); 

    for (size_t j = 0; j < gt_data.size(); j++) {
        double x = gt_data.at(j);
        jacobian.data.at(j).at(0) = x; // dR/d(b1)

        jacobian.data.at(j).at(1) = sign(x) * pow(abs(x), b2) * sin(b3 * pow(abs(x), b4)); // dR/d(b2)

        jacobian.data.at(j).at(2) = b1 * sign(x) * b2 * pow(abs(x), b2) / x
            * sin(b3 * pow(abs(x), b4)); // dR/d(b3)

        jacobian.data.at(j).at(3) = b1 * sign(x) * pow(abs(x), b2) 
            * pow(abs(x), b4) * cos(b3 * pow(abs(x), b4)); // dR/d(b4)

        jacobian.data.at(j).at(4) = b1 * sign(x) * pow(abs(x), b2) 
            * b3 * b4 * pow(abs(x), b4) / x  * cos(b3 * pow(abs(x), b4)); // dR/d(b4) 
    }

    return jacobian;
}
Extended_sin_fit::Extended_sin_fit(int64_t max_steps, double bound, std::vector<double> gt_data, std::vector<double> measurement) {
    double ak = .1;
    double bk = .00001;
    Vector<double> params = Vector<double>(5);
    for (int i = 0; i < params.size; i++) {
        params.data.at(i) = 1;
    }
    this->parameters = params;

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

        Matrix<double> jacobian = this->Jacobian(gt_data);

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
        if (delta_s < bound / 10) {
            break;
        };

        ++steps; 
    }
    printf("residual at end: %f after %d steps \n", S, steps);

}

double Extended_sin_fit::sf_at(double x) {
    return this->parameters.data.at(0) * x + this->parameters.data.at(1) * sign(x) * pow(abs(x), this->parameters.data.at(2)) * sin(this->parameters.data.at(3) * pow(abs(x), this->parameters.data.at(4)));
}

std::vector<double> Extended_sin_fit::sf_calc(std::vector<double> x) {
    std::vector<double> ret;
    for (uint i = 0; i < x.size(); i++) {
        ret.push_back(this->sf_at(x.at(i)));
    }
    return ret;
}

double Extended_sin_fit::sf_at_derivative(double x) {
    double sum = 0;
    double c = this->parameters.data.at(0);
    double b1 = this->parameters.data.at(1);
    double b2 = this->parameters.data.at(2);
    double b3 = this->parameters.data.at(3);
    double b4 = this->parameters.data.at(4);

    sum += c;
    sum += (sign(x) * b1 * b2 * pow(abs(x), b2) * sin(b3 * pow(abs(x), b4))) / x;
    sum += (sign(x) * b1 * pow(abs(x),b2) * b3 * b4 * pow(abs(x),b4) * cos(b3 * pow(abs(x), b4)))/x;
    return sum;
}

double Extended_sin_fit::at(double x, double bound) {
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

std::vector<double> Extended_sin_fit::calc(std::vector<double> x, double bound) {
    std::vector<double> ret;
    for (uint i = 0; i < x.size(); i++) {
        ret.push_back(this->at(x.at(i), bound));
    }
    return ret;
}
