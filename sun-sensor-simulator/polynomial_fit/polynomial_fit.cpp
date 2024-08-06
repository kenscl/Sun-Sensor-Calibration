#include "polynomial_fit.h"
#include "../math/matlib.h"
#include <cstdio>

Matrix<double> Polynomial_Fit::vandermonde(std::vector<double> values){
   Matrix<double> vandi(values.size(), this->degree); 
   for (uint i = 0; i < values.size(); i++) {
       std::vector<double> col;
       for (int j = 0; j < this->degree; j++) {
           col.push_back(pow (values[i], j));
       }
       vandi.data.at(i) = col;
   }
   return vandi;
}

Polynomial_Fit::Polynomial_Fit(int64_t degree, std::vector<double> gt_data, std::vector<double> measurement) {
    this->degree = degree;
    Matrix<double> X(this->degree, this->degree);
    this->weights = std::vector<double>();
    
    X = vandermonde(measurement);
    this->weights = (X.transpose() * X).inverse() * X.transpose() * gt_data;
}

double Polynomial_Fit::at(double x) {
    double sum = 0;
    for (int i = 0; i < this->degree; i++) {
        sum += this->weights.at(i) * pow(x, i);
    }
    return sum;
}

std::vector<double> Polynomial_Fit::calc(std::vector<double> x) {
    std::vector<double> ret;
    for (uint i = 0; i < x.size(); i++) {
        ret.push_back(this->at(x.at(i)));
    }
    return ret;
}
