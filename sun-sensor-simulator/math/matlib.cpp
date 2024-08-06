#include "matlib.h"
#include <cstdio>
#include <stdexcept>

Vector_3D::Vector_3D() {
    this->x = 0;
    this->y = 0;
    this->z = 0;
}

Vector_3D::Vector_3D(double x, double y, double z){
    this->x = x;
    this->y = y;
    this->z = z;
}

Vector_3D::Vector_3D(Vector_3D const &other){
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
}

Vector_3D Vector_3D::cross_product(Vector_3D other) const{
    Vector_3D result;
    result.x = this->y * other.z - this->z * other.y;
    result.y = this->z * other.x - this->x * other.z;
    result.z = this->x * other.y - this->y * other.x;
    return result;
}

double Vector_3D::norm() const{
    double norm;
    norm = sqrt(x*x + y*y + z*z);
    return norm;
}
Vector_3D Vector_3D::normalize() const{
    double norm;
    norm = this->norm();
    if (norm == 0) return Vector_3D(0,0,0);

    Vector_3D result;
    result.x = this->x / norm;
    result.y = this->y / norm;
    result.z = this->z / norm;
    
    return result;
}

double Vector_3D::calculate_angle(Vector_3D other) const{
    double ac;
    Vector_3D v1, v2;
    v1 = this->normalize();
    v2 = other.normalize();
    double res;
    res = v1 * v2;
    ac = acos (res);
    return ac;
}

Matrix_3D Vector_3D::rmat_to(Vector_3D target) const{
    Vector_3D v1 = this->normalize();
    Vector_3D v2 = target.normalize();
    if (v1 == v2) {
        Matrix_3D ret;
        return ret.I();
    }
    Vector_3D vdiff = v2;

    double cos = v1 * v2;
    double length = cos;
    double sin = (v1.cross_product(v2)).norm();


    Matrix_3D rot;
    rot.r[0][0] = cos;
    rot.r[0][1] = - sin;
    rot.r[0][2] = 0;

    rot.r[1][0] = sin;
    rot.r[1][1] = cos;
    rot.r[1][2] = 0;

    rot.r[2][0] = 0; 
    rot.r[2][1] = 0; 
    rot.r[2][2] = 1; 

    vdiff.x -= length * v1.x;
    vdiff.y -= length * v1.y;
    vdiff.z -= length * v1.z;

    vdiff = vdiff.normalize();
    Vector_3D w = v2.cross_product(v1);

    Matrix_3D base_change_inv;

    base_change_inv.r[0][0] = v1.x;
    base_change_inv.r[0][1] = v1.y;
    base_change_inv.r[0][2] = v1.z;

    base_change_inv.r[1][0] = vdiff.x;
    base_change_inv.r[1][1] = vdiff.y;
    base_change_inv.r[1][2] = vdiff.z;

    base_change_inv.r[2][0] = w.x;
    base_change_inv.r[2][1] = w.y;
    base_change_inv.r[2][2] = w.z;

    Matrix_3D base_change;
    base_change = base_change_inv.inverse();
    Matrix_3D rmat = base_change * rot * base_change_inv;

    return rmat;

}

double Vector_3D::distance_to(Vector_3D other) const {
    Vector_3D res;
    res = *this - other;
    return res.norm();
}

void Vector_3D::print() const{
    printf("Vector: \n");
    printf("x: %.5f \n", this->x);
    printf("y: %.5f \n", this->y);
    printf("z: %.5f \n", this->z);
}
Vector_3D Vector_3D::operator*(double value) const{
    Vector_3D result;
    result.x = value * this->x;
    result.y = value * this->y;
    result.z = value * this->z;
    return result;
}

Vector_3D Vector_3D::operator/(double value) const{
    Vector_3D result;
    result.x = this->x / value;
    result.y = this->y / value;
    result.z = this->z / value;
    return result;
}

double Vector_3D::operator*(Vector_3D other) const{
    double res;
    res = this->x * other.x + this->y * other.y + this->z * other.z;
    return res;
}

Vector_3D Vector_3D::operator+(Vector_3D other) const{
    Vector_3D result;
    result.x = this->x + other.x;
    result.y = this->y + other.y;
    result.z = this->z + other.z;
    return result;
}

Vector_3D Vector_3D::operator-(Vector_3D other) const{
    Vector_3D result;
    result.x = this->x - other.x;
    result.y = this->y - other.y;
    result.z = this->z - other.z;
    return result;
}

bool Vector_3D::operator==(Vector_3D other) const{
    if (this->x == other.x && this->y == other.y && this->z == other.z) return true;
    return false;
}

template <typename T>
Vector<T>::Vector(std::vector<T> data) {
    this->data = data;
    this->size = data.size();
}
template <typename T>
Vector<T>::Vector(Vector_3D data) {
    this->data = std::vector<T>();
    this->data.push_back(data.x);
    this->data.push_back(data.y);
    this->data.push_back(data.z);
    this->size = 3;
}

template <typename T>
Vector<T>::Vector(int size) {
   this->data = std::vector<T>(size, 0);
   this->size = size;
}

template <typename T>
int Vector<T>::get_size() const{
    return size;
}
        
template <typename T>
Vector<T> Vector<T>::sub_vector(int start, int end) const {
    if (start < 0 || end > this->size || start > end) {
        throw std::logic_error("Invalid start or end index");
    }
    Vector<T> result(end - start);
    for (int i = start; i < end; i++) {
        result.data.at(i - start) = this->data.at(i);
    }
    return result;
}

template <typename T>
Vector<T> Vector<T>::set_sub_vector(Vector<T> &other, int start, int end) const {
    if (start < 0 || end > this->size || start > end) {
        throw std::logic_error("Invalid start or end index");
    }
    if (end - start != other.size) {
        throw std::logic_error("Subvector size does not match");
    }
    Vector<T> result = *this;
    for (int i = start; i < end; i++) {
        result.data.at(i) = other.data.at(i - start);
    }
    return result;
}

template <typename T>
void Vector<T>::print() const {
    printf("Vector: \n");
    for (int i = 0; i < this->size; i++) {
        printf("%.15f \n", this->data.at(i));
    }
}

template <typename T>
T Vector<T>::norm() const {
    T sum = 0;
    for (int i = 0; i < this->size; ++i) {
        sum += this->data.at(i) * this->data.at(i);
    }
    return sqrt(sum);
}

template <typename T>
Vector<T> Vector<T>::normalize() const {
    if (this->norm() < 1e-10) {
        Vector<T> ret(this->size);
        for (int i = 0; i < this->size; i++) {
            ret.data.at(i) = 0.;
        }
        return ret;
    }
    return (*this)/this->norm();
}

template <typename T>
Vector<T> Vector<T>::operator*(double d) const {
    Vector result(this->size);
    for (int i = 0; i < this->size; i++) {
        result.data.at(i) = this->data.at(i) * d;
    }
    return result;
}

template <typename T>
T Vector<T>::operator*(Vector<T> other) const {
    if (this->size != other.size){
        throw std::logic_error("Vectors not of same length");
    }
    T sum;
    for (int i = 0; i < this->size; i++) {
        sum = sum + this->data.at(i) * other.data.at(i);
    }
    return sum;
}

template <typename T>
Vector<T> Vector<T>::operator/(double d) const {
    Vector result(this->size);
    for (int i = 0; i < this->size; i++) {
        result.data.at(i) = this->data.at(i) / d;
    }
    return result;
}

template <typename T>
Vector<T> Vector<T>::operator+(const Vector<T> other) const {
    if (this->size != other.size){
        throw std::logic_error("Vectors not of same length");
    }
    Vector<T> result(this->size);
    for (int i = 0; i < other.size; i++)  {
        result.data.at(i) = this->data.at(i) + other.data.at(i);
    }
    return result;
}

template <typename T>
Matrix<T> Vector<T>::transpose() const {
    Matrix<T> result(1, this->size);
    for(int i = 0; i < this->size; i++) {
        result.data.at(0).at(i) = this->data.at(i);
    }
    return result;
}
template <typename T>
Vector<T> Vector<T>::operator-(const Vector<T> other) const {
    if (this->size != other.size){
        throw std::logic_error("Vectors not of same length");
    }
    Vector<T> result(this->size);
    for (int i = 0; i < other.size; i++)  {
        result.data.at(i) = this->data.at(i) - other.data.at(i);
    }
    return result;
}

template<typename T>
Vector<T>& Vector<T>::operator=(const Vector<T>& other) {
    if (this != &other) { 
        this->data = other.data;
        this->size = other.size;
    }
    return *this;
}

template<typename T>
Matrix<T> Vector<T>::to_matrix() const {
    Matrix<T> result(this->size, 1);
    for (int i = 0; i < this->size; i++) {
        result.data.at(i).at(0) = this->data.at(i);
    }
    return result;
}
template<typename T>
double Vector<T>::calculate_angle(Vector<T> other) const {
    if (this->size != other.size){
        throw std::logic_error("Vectors not of same length");
    }
    double ac;
    Vector<T> v1, v2;
    v1 = this->normalize();
    v2 = other.normalize();
    double res;
    res = v1 * v2;
    ac = acos (res);
    return ac;
}

template<typename T>
Vector<T> Vector<T>::dot(Vector<T> other) const {
    if (this->size != 3 || other.size != 3) {
        throw std::logic_error("Vectors must be 3 by 3 to calculate cross product");
    }
    Vector<T> result(3);
    result.data.at(0) = this->data.at(1) * other.data.at(2) - this->data.at(2) * other.data.at(1);
    result.data.at(1) = this->data.at(2) * other.data.at(0) - this->data.at(0) * other.data.at(2);
    result.data.at(2) = this->data.at(0) * other.data.at(1) - this->data.at(1) * other.data.at(0);
    return result;
}

template<typename T>
void Vector<T>::householder(Vector<T> &v, double &beta) const {
    if (this->size != v.size) {
        throw std::logic_error("Vectors must be the same size to calculate householder");
    }
    T sigma = 0;
    for (int i = 1; i < this->size; i++) {
        sigma += v.data.at(i) * v.data.at(i);
    }
    T mu = sqrt(v.data.at(0) * v.data.at(0) + sigma);
    if (mu == 0) {
        beta = 0;
        return;
    }
    T v1 = v.data.at(0);
    if (v1 <= 0) {
        v1 = v1 - mu;
    } else {
        v1 = -sigma / (v1 + mu);
    }
    beta = 2 * v1 * v1 / (sigma + v1 * v1);
    v.data.at(0) = 1;
    for (int i = 1; i < this->size; i++) {
        v.data.at(i) = v.data.at(i) / v1;
    }
}
double lerp(double lly, double lry, double llx, double lrx, double x) {
    double distance = abs (lrx - llx);
    if (distance == 0) {
        return (lly + lry) / 2; 
    }
    double weight_left = lly * abs (lrx - x) / distance;
    double weight_right = lry * abs (llx - x) / distance;
    return weight_left + weight_right;
}

double linear_interpolation(double * x, double * y, double value, int len){
    if (x[0] > value) return y[0];
    for (int i = 0; i < len - 1; i++) {
        if (x[i] <= value && x[i+1] >= value) return lerp(y[i], y[i+1], x[i], x[i + 1], value);
    }
    return y[len];
}

double rmse(std::vector<double> u, std::vector<double> v) {
    if (u.size() != v.size()) {
        throw std::logic_error("Datasets must be the same size to calculate rmse");
    }
    double rmse = 0;
    for (int i = 0; i < u.size(); ++i){
        double delta = u.at(i) - v.at(i);
        rmse += delta * delta;
    }
    rmse = sqrt (rmse / u.size());
    return rmse;
}

