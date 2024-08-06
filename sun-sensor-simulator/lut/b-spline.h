#ifndef __B_SPLINE
#define __B_SPLINE
#include <vector>
#include "../math/matlib.h"
class BSpline {
public:
    std::vector<double> x;
    std::vector<double> y;
    Matrix<double> basis_matrix;

public:
    /*
    * This function generates the b-spline of the data that can be used to correct measurements.
    * @param gt_data: The ground truth data that is used to generate the spline.
    * @param measurement: The measurement data that is used to generate the spline.
    * @param num_points: The number of knots
    */
    BSpline(std::vector<double> gt_data, std::vector<double> measurement, int num_points);

    /*
     * This function generates the coefficients of the b-spline of the data.
     * @param points: The points that are used to generate the spline.
     */
    std::vector<std::vector<double>> generate_coefficients(std::vector<double> x, std::vector<double> y);

    /*
    * This function evalutes the spline at a given point.
    * @param x: The point at which the spline is evaluated.
    */
    std::vector<double> get_value(double x);
    std::vector<std::vector<double>> get_value(std::vector<double> x);

    double at(double x);
    std::vector<double> calc(std::vector<double> x);
    
    double basis_function(int i, int p, double t, std::vector<double> knots);


};
#endif // __B_SPLINE
