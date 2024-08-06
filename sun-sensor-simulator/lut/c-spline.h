#ifndef __C_SPLINE
#define __C_SPLINE
#include <vector>
#include "../math/matlib.h"
class CSpline {
public:
    std::vector<std::vector<double>> coefficients; 
    std::vector<double> y;

public:
    /*
    * This function generates the natural cubic spline interpolation of the data that can be used to correct measurements.
    * @param gt_data: The ground truth data that is used to generate the spline.
    * @param measurement: The measurement data that is used to generate the spline.
    * @param num_points: The number of knots
    */
    CSpline(std::vector<double> gt_data, std::vector<double> measurement, int num_points);

    /*
     * This function generates the coefficients of the natural cubic spline interpolation of the data.
     * @param points: The points that are used to generate the spline.
     */
    std::vector<std::vector<double>> generate_coefficients(std::vector<double> x, std::vector<double> y);

    /*
    * This function evalutes the spline at a given point.
    * @param x: The point at which the spline is evaluated.
    */
    double get_value(double x);
    std::vector<double> get_value(std::vector<double> x);

    double at(double x);
    std::vector<double> calc(std::vector<double> x);

};
#endif // __C_SPLINE
