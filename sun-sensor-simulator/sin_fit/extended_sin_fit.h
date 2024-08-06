#ifndef __Extended_SIN_FIT
#define __Extended_SIN_FIT
#include <cstdint>
#include <../math/matlib.h>
class Extended_sin_fit
{
    public:
        Vector<double> parameters;

    private:
        Matrix<double> Jacobian(const std::vector<double> &gt_data);

    public:
        /*
         * Constructor, generates the sin fit.
         * gt_data is the actual data
         * measurement is the measurement
         * the function that the data is fitted to is f(x) = c * x + b:1 * sign(x) * abs (x)^b_2 * sin (b_3 * abs(x)^b_4)
         * see: https://de.wikipedia.org/wiki/Levenberg-Marquardt-Algorithmus, minimisaton happens at the origin
         *
         * NOTE: This function does not correct the measurements but just generate a function of form (1),
         *       that maps the gt data to the measurements.
         */
        Extended_sin_fit(int64_t max_steps, double bound, std::vector<double> gt_data, std::vector<double> measurement);
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

