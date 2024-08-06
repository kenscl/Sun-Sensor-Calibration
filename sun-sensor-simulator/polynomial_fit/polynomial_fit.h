#ifndef __POL_FIT
#define __POL_FIT
#include "../math/matlib.h"
#include <vector>

class Polynomial_Fit {
    public:
        int64_t degree;
        std::vector<double> weights;
    private: 
       Matrix<double>  vandermonde(std::vector<double> values);
    public:
        /*
         * Constructor, generates the polynomial fit.
         * gt_data is the actual data
         * measurement is the measurement
         * degree is the degree of the polynomial, so degree = 2 would be a1 + a2 * x + a3 * x^2 
         */
        Polynomial_Fit(int64_t degree, std::vector<double> gt_data, std::vector<double> measurement);
        /*
         * evaluate the polynomial at a certain point
         */
        double at(double x);
        /*
         * returns the values in x evaluated by at(double x);
         */
        std::vector<double> calc(std::vector<double> x);
};
#endif
