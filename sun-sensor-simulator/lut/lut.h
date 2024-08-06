#ifndef __LUT
#define __LUT
#include "../math/matlib.h"
#include <vector>

class LUT {
    public:
        std::vector<double> parameters;
        std::vector<double> measurement;

    public:
        /*
         * Constructor, generates the Look up table.
         * gt_data is the actual data
         * measurement is the measurement
         *
         * This function generates a look up table from the data points and then
         * returns the corrected data through linear interpolation
         */
        LUT(std::vector<double> gt_data, std::vector<double> measurement, int num_points);
        /*
         * evaluate the lut at x throuth linear interpolation
         */
        double at(double x);
        /*
         * returns the values in x evaluated by at(double x);
         */
        std::vector<double> calc(std::vector<double> x);
};
#endif
