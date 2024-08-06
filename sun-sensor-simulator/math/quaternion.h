#ifndef __QUATERNION
#define __QUATERNION
#include "math.h"
#include "matrix.h"
#include <cmath>

template <typename T>
class Vector;

template <typename T>
class Matrix;

template <typename T>
class Quaternion {
    public:
        T q, i, j, k;
    public: 
        Quaternion();
        Quaternion(T q, T i, T j, T k);
        Quaternion(T angle, const Vector<T>& axis);
        Quaternion(const Matrix<T>& rotation_matrix);
        Quaternion(const Vector<T>& euler_angles);

        Quaternion<T> operator*(const Quaternion<T>& other) const;
        Quaternion<T> operator+(const Quaternion<T>& other) const;
        Quaternion<T> operator*(const T& nbr) const;

        Quaternion<T> operator*(const Matrix<T>& other) const;

        Quaternion<T>& operator=(const Quaternion<T>& other);
        Quaternion<T> conjugate() const;
        Vector<T> to_rpy() const;
        Matrix<T> to_rotation_matrix() const;

        T norm() const;
        Quaternion<T> normalize() const;

        void print() const;
        Quaternion<T> get_diff(const Quaternion<T>& other) const;
};

template class Quaternion<double>;
#endif
