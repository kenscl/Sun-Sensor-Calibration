#ifndef __SIN_FIT
#define __SIN_FIT
#include "../math/matlib.h"
#include <vector>

class Sin_Fit {
    public:
        int64_t degree;
        Vector<double> parameters;

    private:
        Matrix<double> Jacobian(int64_t degree, const std::vector<double> &gt_data);

    public:
        /*
         * Constructor, generates the sin fit.
         * gt_data is the actual data
         * measurement is the measurement
         * degree is the number of terms in the sine, so degree = 2 would be b_1 * sin (g_1 * x) + b_2 * sin (g_2 * x)  (1)
         * see: https://de.wikipedia.org/wiki/Levenberg-Marquardt-Algorithmus, minimisaton happens at the origin
         *
         * NOTE: This function does not correct the measurements but just generate a function of form (1),
         *       that maps the gt data to the measurements.
         */
        Sin_Fit(int64_t degree, int64_t max_steps, double bound, std::vector<double> gt_data, std::vector<double> measurement);
        /*
         * evaluate the (1) at a certain point
         */
        double sf_at(double x);
        /*
         * returns the values in x evaluated by at(double x);
         */
        double sf_at_derivative(double x);
        /*
         * returns the derivative;
         */
        std::vector<double> sf_calc(std::vector<double> x);
        /*
         * evaluate the inverted form of (1) at a certain point
         */
        double at(double x, double bound);
        /*
         * returns the values in x evaluated by at(double x);
         */
        std::vector<double> calc(std::vector<double> x, double bound);
};
#endif
