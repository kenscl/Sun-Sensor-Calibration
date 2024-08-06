#include "quaternion.h"
#include "matlib.h"
#include <system_error>

template <typename T>
class Vector;

template <typename T>
Quaternion<T>::Quaternion() {
    this->q = 0;
    this->i = 0;
    this->j = 0;
    this->k = 0;
}

template <typename T>
Quaternion<T>::Quaternion(T q, T i, T j, T k) {
    this->q = q;
    this->i = i;
    this->j = j;
    this->k = k;
}

template <typename T>
Quaternion<T>::Quaternion(T angle, const Vector<T>& axis) {
    this->q = cos(angle / 2);
    this->i = axis.data.at(0) * sin(angle / 2);
    this->j = axis.data.at(1) * sin(angle / 2);
    this->k = axis.data.at(2) * sin(angle / 2);
}

template <typename T>
Quaternion<T>::Quaternion(const Matrix<T>& rotation_matrix) {
    if (rotation_matrix.rows != 3 || rotation_matrix.cols != 3) {
        throw std::invalid_argument("rotation matrix must be 3x3");
    }
    Quaternion<T> ret = Quaternion<T>(1, 0, 0, 0);

    T trace = rotation_matrix.data.at(0).at(0) + rotation_matrix.data.at(1).at(1) + rotation_matrix.data.at(2).at(2);

    if (trace > 0) {
        T s = 0.5 / sqrt(trace + 1);
        this->q = 0.25 / s;
        this->i = (rotation_matrix.data.at(2).at(1) - rotation_matrix.data.at(1).at(2)) * s;
        this->j = (rotation_matrix.data.at(0).at(2) - rotation_matrix.data.at(2).at(0)) * s;
        this->k = (rotation_matrix.data.at(1).at(0) - rotation_matrix.data.at(0).at(1)) * s;
    } else {
        if (rotation_matrix.data.at(0).at(0) > rotation_matrix.data.at(1).at(1) && rotation_matrix.data.at(0).at(0) > rotation_matrix.data.at(2).at(2)) {
            T s = 2 * sqrt(1 + rotation_matrix.data.at(0).at(0) - rotation_matrix.data.at(1).at(1) - rotation_matrix.data.at(2).at(2));
            this->q = (rotation_matrix.data.at(2).at(1) - rotation_matrix.data.at(1).at(2)) / s;
            this->i = 0.25 * s;
            this->j = (rotation_matrix.data.at(0).at(1) + rotation_matrix.data.at(1).at(0)) / s;
            this->k = (rotation_matrix.data.at(0).at(2) + rotation_matrix.data.at(2).at(0)) / s;
        } else if (rotation_matrix.data.at(1).at(1) > rotation_matrix.data.at(2).at(2)) {
            T s = 2 * sqrt(1 + rotation_matrix.data.at(1).at(1) - rotation_matrix.data.at(0).at(0) - rotation_matrix.data.at(2).at(2));
            this->q = (rotation_matrix.data.at(0).at(2) - rotation_matrix.data.at(2).at(0)) / s;
            this->i = (rotation_matrix.data.at(0).at(1) + rotation_matrix.data.at(1).at(0)) / s;
            this->j = 0.25 * s;
            this->k = (rotation_matrix.data.at(1).at(2) + rotation_matrix.data.at(2).at(1)) / s;
        } else {
            T s = 2 * sqrt(1 + rotation_matrix.data.at(2).at(2) - rotation_matrix.data.at(0).at(0) - rotation_matrix.data.at(1).at(1));
            this->q = (rotation_matrix.data.at(1).at(0) - rotation_matrix.data.at(0).at(1)) / s;
            this->i = (rotation_matrix.data.at(0).at(2) + rotation_matrix.data.at(2).at(0)) / s;
            this->j = (rotation_matrix.data.at(1).at(2) + rotation_matrix.data.at(2).at(1)) / s;
            this->k = 0.25 * s;
        }
    }
}
template <typename T>
Quaternion<T>::Quaternion(const Vector<T>& euler_angles) {
    this->q = 0;
    this->i = euler_angles.data.at(0);
    this->j = euler_angles.data.at(1);
    this->k = euler_angles.data.at(2);
}

template <typename T>
Quaternion<T> Quaternion<T>::operator*(const Quaternion<T>& other) const {
    return Quaternion<T>(
        q * other.q - i * other.i - j * other.j - k * other.k,  
        q * other.i + i * other.q + j * other.k - k * other.j,  
        q * other.j - i * other.k + j * other.q + k * other.i,  
        q * other.k + i * other.j - j * other.i + k * other.q   
    );
}

template <typename T>
Quaternion<T> Quaternion<T>::operator+(const Quaternion<T>& other) const {
    return Quaternion<T>(
        q + other.q,
        i + other.i,
        j + other.j,
        k + other.k
    );
}

template <typename T>
Quaternion<T> Quaternion<T>::operator*(const Matrix<T>& other) const {
    Vector<T> v = Vector<T>(4);
    v.data.at(0) = this->q;
    v.data.at(1) = this->i;
    v.data.at(2) = this->j;
    v.data.at(3) = this->k;
    v = other * v;
    return Quaternion<T>(v.data.at(0), v.data.at(1), v.data.at(2), v.data.at(3));
}

template <typename T>
Quaternion<T> Quaternion<T>::operator*(const T& nbr) const {
    return Quaternion<T>(q * nbr, i * nbr, j * nbr, k * nbr);
}
template <typename T>
Quaternion<T> Quaternion<T>::conjugate() const {
    return Quaternion<T>(q, -i, -j, -k);
}

template <typename T>
Vector<T> Quaternion<T>::to_rpy() const {
    Vector<T> result(3);
    result.data.at(0) = atan2(2 * (q * i + j * k), 1 - 2 * (i * i + j * j));
    result.data.at(1) = asin(2 * (q * j - k * i));
    result.data.at(2) = atan2(2 * (q * k + i * j), 1 - 2 * (j * j + k * k));
    return result;
}

template <typename T>
Matrix<T> Quaternion<T>::to_rotation_matrix() const {
    Matrix<T> result(3, 3);
    result.data.at(0).at(0) = 1 - 2 * (j * j + k * k);
    result.data.at(0).at(1) = 2 * (i * j - k * q);
    result.data.at(0).at(2) = 2 * (i * k + j * q);
    result.data.at(1).at(0) = 2 * (i * j + k * q);
    result.data.at(1).at(1) = 1 - 2 * (i * i + k * k);
    result.data.at(1).at(2) = 2 * (j * k - i * q);
    result.data.at(2).at(0) = 2 * (i * k - j * q);
    result.data.at(2).at(1) = 2 * (j * k + i * q);
    result.data.at(2).at(2) = 1 - 2 * (i * i + j * j);
    return result;
}

template <typename T>
Quaternion<T>& Quaternion<T>::operator=(const Quaternion<T>& other) {
    if (this != &other) {  
        this->q = other.q;
        this->i = other.i;
        this->j = other.j;
        this->k = other.k;
    }
    return *this;
}


template <typename T>
T Quaternion<T>::norm() const {
    return sqrt(q * q + i * i + j * j + k * k);
}
template <typename T>
Quaternion<T> Quaternion<T>::normalize() const {
    T n = this->norm();
    if (n == 0) {
        return Quaternion<T>(0, 0, 0, 0);
    }
    return Quaternion<T>(q / n, i / n, j / n, k / n);
}

template <typename T>
void Quaternion<T>::print() const {
    printf("quaternion: \n");
    printf("q %.5f \n", this->q);
    printf("i %.5f \n", this->i);
    printf("j %.5f \n", this->j);
    printf("k %.5f \n", this->k);
}

template <typename T>
Quaternion<T> Quaternion<T>::get_diff(const Quaternion<T>& other) const {
    double norm = other.norm();
    if (norm == 0) {
        return Quaternion<T>(*this);
    }
    norm *= norm;
    norm = 1 / norm;
    Quaternion<T> q = other.conjugate() * norm;
    return *this * q;
}
