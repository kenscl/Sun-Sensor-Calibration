#ifndef __MATLIB
#define __MATLIB
#include "math.h"
#include <vector>

// Math includes:
#include "matrix.h"
#include "quaternion.h"

// constants
# define PI                 3.14159265359
# define C                  299792458.0
# define PLANCK_CONSTANT    6.62607015e-34
# define BOLTZMANN_CONSTANT 1.380649e-23

#define R2D     180 / M_PI
#define D2R     M_PI / 180

class Vector_3D;
class Matrix_3D;
template <typename T>
class Vector;

class Vector_3D {
    public: 
        double x,y,z;

    public: 
        Vector_3D();
        Vector_3D(double x, double y, double z);
        Vector_3D(Vector_3D const &other);

        Vector_3D cross_product(Vector_3D other) const;
        double norm() const;
        Vector_3D normalize() const;
        void print() const;
        double calculate_angle(Vector_3D other) const;
        Matrix_3D rmat_to(Vector_3D target) const;
        double distance_to(Vector_3D other) const;

        Vector_3D operator*(double value) const;
        Vector_3D operator/(double value) const;
        bool operator==(Vector_3D other) const;

        double operator*(Vector_3D other) const;
        Vector_3D operator+(Vector_3D other) const;
        Vector_3D operator-(Vector_3D other) const;
};

template <typename T>
class Vector {
    public:
        std::vector<T> data;
        int size;
    public:
        Vector(std::vector<T> data);
        Vector(Vector_3D data);
        Vector(int size);
        Vector() : size(0) {}

        int get_size() const;
        double calculate_angle(Vector<T> other) const;
        void print() const;
        T norm() const;
        Vector<T> normalize() const;
        Vector<T> sub_vector(int start, int end) const;
        Vector<T> set_sub_vector(Vector<T> &other, int start, int end) const;
        Vector<T> dot(Vector<T> other) const;
        void householder(Vector<T> &v, double &beta) const;
        Matrix<T> transpose() const;
        Matrix<T> to_matrix() const;
        
        Vector<T> operator*(double d) const;
        T operator*(Vector<T> other) const;
        Vector<T> operator/(double d) const;
        Vector<T> operator-(const Vector<T> other) const;
        Vector<T> operator+(const Vector<T> other) const;
        Vector<T>& operator=(const Vector<T>& other);
};

double lerp(double lly, double lry, double llx, double lrx, double x);
double linear_interpolation(double * x, double * y, double value, int len);
double rmse(std::vector<double> u, std::vector<double> v);

inline double sign(double x) {
    if (x < 0) return -1;
    return 1;
}

inline void convert_to_deg(std::vector<double> &data) {
    for (int i = 0; i < data.size(); i++) {
        data.at(i) = data.at(i) * R2D;
    }
}

inline double rmse_45(std::vector<double> data, std::vector<double> gt) {
    double sum = 0;
    for (int i = 0; i < data.size(); i++) {
        if (gt.at(i) < -45 * D2R || gt.at(i) > 45 * D2R) {
            continue;
        }
        sum += pow(data.at(i) - gt.at(i), 2);
    }
    return sqrt(sum / data.size());
}

inline std::vector<double> select_evenly_spaced_points(const std::vector<double>& points, int n) {
    std::vector<double> selectedPoints;
    int totalPoints = points.size();
    if (n <= 0 || n > totalPoints) {
        return points;
    }
    
    double step = static_cast<double>(totalPoints - 1) / (n - 1);
    
    for (int i = 0; i < n; ++i) {
        int index = static_cast<int>(i * step);
        selectedPoints.push_back(points[index]);
    }
    
    return selectedPoints;
}

inline void print_vector(std::vector<double> vec) {
    for (uint i = 0; i < vec.size(); ++i) {
        printf("%f\n", vec.at(i));
    }
}

inline std::vector<double> get_derivative(std::vector<double> x, std::vector<double> y) {
    std::vector<double> ret;
    for (int i = 0; i < x.size() - 1; i++) {
        ret.push_back((y.at(i + 1) - y.at(i)) / (x.at(i + 1) - x.at(i)));
    }
        ret.push_back(ret.back());
    return ret;
}

template <typename  TYPE>
inline int concatenate_axis_1(Matrix<TYPE> &A, Matrix<TYPE> &B, Matrix<TYPE> & return_value){
    /*
     * return matrix is given in return_value
     * If error occurs 1 is returned, else 0
     */
    if (A.rows != B.rows) {
        return 1;
    }
    for (int i = 0; i < (int) A.rows; i++){
        for (int j = 0; j < (int) A.cols; j++) {
            return_value.data[i][j] = A.data[i][j];
        }
    }
    for (int i = 0; i < (int) A.rows; i++){
        for (int j = (int) A.cols; j < (int) A.cols + (int) B.cols; j++) {
            return_value.data[i][j] = B.data[i][j - (int) A.cols];
        }
    }
    return 0;
}
template <typename  TYPE>
inline int concatenate_axis_0(Matrix<TYPE> &A, Matrix<TYPE> &B, Matrix<TYPE>& return_value){
    /*
     * return matrix is given in return_value
     * If error occurs 1 is returned, else 0
     */

    if (A.cols != B.cols) {
        return  1;
    }
    for (int i = 0; i < (int) A.rows; i++){
        for (int j = 0; j < (int) A.cols; j++) {
            return_value.data[i][j] = A.data[i][j];
        }
    }
    for (int i = (int) A.rows; i < (int) A.rows + (int) B.rows; i++){
        for (int j = 0; j < (int) A.cols; j++) {
            return_value.data[i][j] = B.data[i - (int) A.rows][j];
        }
    }
    return 0;

}

template class Vector<double>;
#endif

