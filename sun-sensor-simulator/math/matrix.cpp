#include "matrix.h"
#include "matlib.h"
#include <cmath>
#include <vector>

Matrix_3D::Matrix_3D(){
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            this->r[i][j] = 0;
        }
    }
}


Matrix_3D Matrix_3D::I() const{
    Matrix_3D ret;
    ret.r[0][0] = 1;
    ret.r[0][1] = 0;
    ret.r[0][2] = 0;

    ret.r[1][0] = 0;
    ret.r[1][1] = 1;
    ret.r[1][2] = 0;

    ret.r[2][0] = 0;
    ret.r[2][1] = 0;
    ret.r[2][2] = 1;
    return ret;
}

Vector_3D Matrix_3D::operator*(Vector_3D v) const{
    Vector_3D result;
    result.x = this->r[0][0] * v.x + this->r[0][1] * v.y + this->r[0][2] * v.z;
    result.y = this->r[1][0] * v.x + this->r[1][1] * v.y + this->r[1][2] * v.z;
    result.z = this->r[2][0] * v.x + this->r[2][1] * v.y + this->r[2][2] * v.z;
    return result;
}

Matrix_3D& Matrix_3D::operator=(const Matrix_3D& other){
        if (this != &other) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    r[i][j] = other.r[i][j];
                }
            }
        }
        return *this;
}
Matrix_3D Matrix_3D::operator*(const Matrix_3D& other) const{
    Matrix_3D result;
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                result.r[i][j] += this->r[i][k] * other.r[k][j];
            }
        }
    }
    return result;
}

Matrix_3D Matrix_3D::operator*(const double scalar) const {
    Matrix_3D result;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result.r[i][j] = r[i][j] * scalar;
        }
    }
    return result;
}
Matrix_3D Matrix_3D::operator+(const Matrix_3D& other) const {
        Matrix_3D result;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                result.r[i][j] = r[i][j] + other.r[i][j];
            }
        }
        return result;
    }

void Matrix_3D::print() const{
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%f ", r[i][j]);
            if (j == 2) printf("\n");
        }
    }
}

double Matrix_3D::determinant() const{
    double det = this->r[0][0] * (this->r[1][1] * this->r[2][2] - this->r[1][2] * this->r[2][1]) -
        this->r[0][1] * (this->r[1][0] * this->r[2][2] - this->r[1][2] * this->r[2][0]) +
        this->r[0][2] * (this->r[1][0] * this->r[2][1] - this->r[1][1] * this->r[2][0]);
    return det;
}
Matrix_3D Matrix_3D::inverse() const{
        double det = this->determinant();

        if (det == 0.0) {
            printf("matrix is singular no inversion possible \n");
        }

        Matrix_3D inv;

        inv.r[0][0] = (this->r[1][1] * this->r[2][2] - this->r[1][2] * this->r[2][1]) / det;
        inv.r[0][1] = (this->r[0][2] * this->r[2][1] - this->r[0][1] * this->r[2][2]) / det;
        inv.r[0][2] = (this->r[0][1] * this->r[1][2] - this->r[0][2] * this->r[1][1]) / det;

        inv.r[1][0] = (this->r[1][2] * this->r[2][0] - this->r[1][0] * this->r[2][2]) / det;
        inv.r[1][1] = (this->r[0][0] * this->r[2][2] - this->r[0][2] * this->r[2][0]) / det;
        inv.r[1][2] = (this->r[0][2] * this->r[1][0] - this->r[0][0] * this->r[1][2]) / det;

        inv.r[2][0] = (this->r[1][0] * this->r[2][1] - this->r[1][1] * this->r[2][0]) / det;
        inv.r[2][1] = (this->r[0][1] * this->r[2][0] - this->r[0][0] * this->r[2][1]) / det;
        inv.r[2][2] = (this->r[0][0] * this->r[1][1] - this->r[0][1] * this->r[1][0]) / det;

        return inv;
}


/*
 * General Matrix Class
 */

template <typename T>
Matrix<T>::Matrix() {}
template <typename T>
Matrix<T>::Matrix(int n, int m) : rows(n), cols(m), data(n, std::vector<T>(m, 0)) {}

template <typename T>
Matrix<T> Matrix<T>::Identity(int size) {
    Matrix<T> ret(size, size);
    for (int i = 0; i < size; ++i) {
        ret.data[i][i] = 1;
    }
    return ret;
}

template <typename T>
std::vector<T> Matrix<T>::operator*(const std::vector<T>& vec) const {
    if (vec.size() != cols) {
        printf("Vector size must match the number of matrix columns \n");
    }
    std::vector<T> result(rows, 0);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i] += data[i][j] * vec[j];
        }
    }
    return result;
}

template <typename T>
Vector<T> Matrix<T>::operator*(const Vector<T>& vec) const {
    if (vec.get_size() != cols) {
        throw std::logic_error("Vector size must match the number of matrix columns \n");
    }
    Vector<T> result(rows);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.data.at(i) += data.at(i).at(j) * vec.data.at(j);
        }
    }
    return result;
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other) {
    if (this != &other) {
        rows = other.rows;
        cols = other.cols;
        data = other.data;
    }
    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& other) const {
    if (cols != other.rows) {
        throw std::logic_error("Matrix dimensions must agree for multiplication\n");
    }
    Matrix<T> result(rows, other.cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < other.cols; ++j) {
            for (int k = 0; k < cols; ++k) {
                result.data[i][j] += data[i][k] * other.data[k][j];
            }
        }
    }
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const T scalar) const {
    Matrix<T> result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.data[i][j] = data[i][j] * scalar;
        }
    }
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) const {
    if (rows != other.rows || cols != other.cols) {
        printf("dimensions: %dx%d %dx%d \n", rows, cols, rows, other.rows);
        throw std::logic_error("Matrix dimensions must agree for addition");
    }
    Matrix<T> result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.data[i][j] = data[i][j] + other.data[i][j];
        }
    }
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& other) const {
    if (rows != other.rows || cols != other.cols) {
        throw std::logic_error("Matrix dimensions must agree for addition");
    }
    Matrix<T> result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.data[i][j] = data[i][j] - other.data[i][j];
        }
    }
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::I() {
    if (rows != cols) {
        throw std::logic_error("Identity matrix must be square");
    }
    Matrix<T> result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        result.data[i][i] = 1;
    }
    return result;
}

template <typename T>
void Matrix<T>::print() const {
    printf("matrix: \n");
    for (const auto& row : data) {
        for (const auto& elem : row) {
            printf("%.5f ", elem) ;
        }
        printf("\n");
    }
}

template <typename T>
Matrix<T> Matrix<T>::submatrix(int excludeRow, int excludeCol) const {
    if (excludeRow < 0 || excludeRow >= rows || excludeCol < 0 || excludeCol >= cols) {
        throw std::out_of_range("Row or column to exclude is out of range");
    }
    Matrix<T> sub(rows - 1, cols - 1);
    int subi = 0;
    for (int i = 0; i < rows; ++i) {
        if (i == excludeRow) continue;
        int subj = 0;
        for (int j = 0; j < cols; ++j) {
            if (j == excludeCol) continue;
            sub.data[subi][subj] = data[i][j];
            ++subj;
        }
        ++subi;
    }
    return sub;
}

template <typename T>
T Matrix<T>::determinant() const {
    if (rows != cols) {
        throw std::logic_error("Determinant can only be calculated for square matrices");
    }
   int n = rows;
    std::vector<std::vector<T>> temp = data; 

    T det = 1;
    for (int i = 0; i < n; ++i) {
        T maxElem = temp[i][i];
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(temp[k][i]) > std::abs(maxElem)) {
                maxElem = temp[k][i];
                maxRow = k;
            }
        }

        if (maxRow != i) {
            std::swap(temp[i], temp[maxRow]);
            det *= -1; 
        }
        if (temp[i][i] == 0) {
            return 0;
        }

        for (int k = i + 1; k < n; ++k) {
            T factor = temp[k][i] / temp[i][i];
            for (int j = i; j < n; ++j) {
                temp[k][j] -= factor * temp[i][j];
            }
        }
        det *= temp[i][i];
    }
    return det;

}

template <typename T>
Matrix<T> Matrix<T>::adjugate() const {
    if (rows != cols) {
        throw std::logic_error("Adjugate can only be calculated for square matrices");
    }
    Matrix<T> adj(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            adj.data[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * submatrix(i, j).determinant();
        }
    }
    return adj;
}

template <typename T>
Matrix<T> Matrix<T>::inverse() const {
if (rows != cols) {
        throw std::logic_error("Inverse can only be calculated for square matrices");
    }
    T det = determinant();
    if (det == 0) {
        throw std::logic_error("Matrix is singular and cannot be inverted.");
    }
    Matrix<T> adj = adjugate();
    return adj * (1 / det);
}

template <typename T>
Matrix<T> Matrix<T>::transpose() const {
    Matrix<T> result(cols, rows);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.data[j][i] = data[i][j];
        }
    }
    return result;
}

template <typename T>
Vector<T> Matrix<T>::get_col(int col) const {
    if (col < 0 || col >= cols) {
        throw std::out_of_range("Column index out of range");
    }

    Vector<T> column(rows);
    for (int i = 0; i < rows; ++i) {
        column.data.at(i) = this->data.at(i).at(col);
    }
    return column;
}

template <typename T>
void Matrix<T>::swap_rows(int row1, int row2) {
    if (row1 < 0 || row1 >= rows || row2 < 0 || row2 >= rows) {
        throw std::out_of_range("Row index out of range");
    }

    for (int j = 0; j < cols; ++j) {
        std::swap(this->data.at(row1).at(j), this->data.at(row2).at(j));
    }
}

template <typename T>
void Matrix<T>::swap_cols(int col1, int col2) {
    if (col1 < 0 || col1 >= cols || col2 < 0 || col2 >= cols) {
        throw std::out_of_range("Column index out of range");
    }

    for (int i = 0; i < rows; ++i) {
        std::swap(this->data.at(i).at(col1), this->data.at(i).at(col2));
    }
}

template <typename T>
Vector<T> Matrix<T>::diag() const {
    if (rows != cols) {
        throw std::logic_error("Diagonal can only be calculated for square matrices");
    }
    Vector<T> diagonal(rows);
    for (int i = 0; i < rows; ++i) {
        diagonal.data.at(i) = data.at(i).at(i);
    }
    return diagonal;
}

template <typename T>
void Matrix<T>::diag(Vector<T> v){
    if (rows != cols) {
        throw std::logic_error("Diagonal can only be calculated for square matrices");
    }
    for (int i = 0; i < rows; ++i) {
        this->data[i][i] = v.data[i];
    }
}


template <typename T>
void Matrix<T>::set_col(int col, Vector<T> v) {
    if (col < 0 || col >= cols) {
        throw std::out_of_range("Column index out of range");
    }
    if (v.get_size() != rows) {
        throw std::logic_error("Vector size must match the number of matrix rows");
    }

    for (int i = 0; i < rows; ++i) {
        data.at(i).at(col) = v.data.at(i);
    }
}


template <typename T>
Vector<T> gaussian_elimination(Matrix<T> A, Vector<T> b) {
       int n = A.size();

       for (int i = 0; i < n - 1; ++i) {
           int pivot_row = i;
           for (int k = i + 1; k < n; ++k) {
               if (std::abs(A.data[k][i]) > std::abs(A.data[pivot_row][i])) {
                   pivot_row = k;
               }
           }
           if (pivot_row != i) {
               A = A.swap_rows(i, pivot_row);
               b.data.swap(i, pivot_row);
           }

           for (int k = i + 1; k < n; ++k) {
               double factor = A.data[k][i] / A.data[i][i];
               for (int j = i; j < n; ++j) {
                   A.data[k][j] -= factor * A.data[i][j];
               }
               b.data[k] -= factor * b.data[i];
           }
       }

       // Back substitution
       Vector<T> x(n);
       for (int i = n - 1; i >= 0; --i) {
           double sum = 0.0;
           for (int j = i + 1; j < n; ++j) {
               sum += A.data[i][j] * x.data[j];
           }
           x.data[i] = (b.data[i] - sum) / A.data[i][i];
       }

       return x;

}

template <typename T>
void Matrix<T>::QR_decomp(Matrix<T> &A, Matrix<T> &Q, Matrix<T> &R) {
    if (A.cols != A.rows) {
        throw std::logic_error("A is not suited to QR decomposition");
    }
    int n = A.rows;
    Matrix<double> A_cp = A;
    Q = Matrix<double>::Identity(n);

    for (int k = 0; k < n; ++k) {
        // Find Householder reflector
        Vector<double> x = A_cp.get_col(k);
        for (int i = 0; i < k; ++i) {
            x.data[i] = 0; // Zero out the upper part
        }
        double alpha = x.norm();
        if (A.data[k][k] < 0) {
            alpha = abs(alpha);
        }
        else {
            alpha = -abs(alpha);
        }

        Vector<double> e(n);
        e.data[k] = 1.0;
        Vector<double> v = (x + e * alpha).normalize();

        Matrix<double> v_mat(n, 1);
        for (int i = 0; i < n; ++i) {
            v_mat.data[i][0] = v.data[i];
        }
        Matrix<double> H = Matrix<double>::Identity(n) - (v_mat * v_mat.transpose()) * 2.0;

        A_cp = H * A_cp;
        Q = Q * H.transpose();
    }

    R = A_cp;
    Q = Q; // Q should be orthogonal
}

template <typename T>
void Matrix<T>::QR_decomp_symmetric(Matrix<T> &A_, Matrix<T> &Q, Matrix<T> &R) {
    Matrix<T> A = A_;
    int n = A.rows;
    for (int i = 0; i < n - 2; ++i) {
        double beta = 0;
        Vector<T> to_household(n - i - 1);
        for (int j = i + 1; j < n; ++j) {
            to_household.data[j-i-1] = A.data[i][j];
        }
        Vector<T> v(to_household.size);

        to_household.householder(v, beta);
        Matrix<T> A_part(n-i-1, n-i-1);
        Vector<T> v_part(n-i-1);
        for (int j = i + 1; j < A.rows; ++j) {
            for (int k = i + 1; k < A.cols; ++k) {
                A_part.data[j - i - 1][k - i - 1] = A.data[j][k];
                v_part.data[j - i - 1] = A.data[j][i];
            }
        }

        Vector<T> p =  A_part * v * beta;
        printf("there \n");
        Vector<T> w = v - (p.transpose() * v * beta * 0.5) * v;
        A.data[i+1][i] = v_part.norm();
        A.data[i][i+1] = A.data[i+1][i];
        Matrix<T> A_new(n-i-1, n-i-1);
        for (int j = i + 1; j < A.rows; ++j) {
            for (int k = i + 1; k < A.cols; ++k) {
                A_new.data[j - i - 1][k - i - 1] = A.data[j][k];
            }
        }
        A_new = A_new - v.to_matrix() * w.transpose() - w.to_matrix() * v.transpose();
        for (int j = i + 1; j < A.rows; ++j) {
            for (int k = i + 1; k < A.cols; ++k) {
                A.data[j][k] = A_new.data[j - i - 1][k - i - 1];
            }
        }
    }

}

template <typename T>
void Matrix<T>::QR_column_pivoting(Matrix<T> &input, Matrix<T> &Q, Matrix<T> &R) {
    Matrix<T> A = input;
    std::vector<T> c;
    for (int i = 0; i < A.cols; ++i) {
        Vector<T> col(A.rows);
        for (int j = 0; j < A.rows; ++j) {
            col.data[j] = A.data[j][i];
        }
        c.push_back((A.transpose() * A).data[0][0]);
    }
    int r = 0;
    double t = - MAXFLOAT;
    for (int i = 0; i < c.size(); ++i) {
        if (c[i] > t) {
            t = c[i];
        }
    }

    while (t > 0 && r < A.cols) {
        r += 1;
        int k = 0;
        for (int i = r; i < c.size(); ++i) {
            if (c[i] == t) {
                k = i;
                break;
            }
        }
        A.swap_cols(r, k);
        T temp_r = c[r];
        c[r] = c[k];
        c[k] = temp_r;

        Vector<T> to_household;
        for (int i = r + 1; i < A.cols; ++i) {
            to_household.data.push_back(A.data[r][i]);
        }
        double beta = 0;
        Vector<T> v;
        to_household.householder(v, beta);
        Matrix<T> H = Matrix<T>::Identity(A.cols) - v.to_matrix() * v.transpose() * beta;
        Matrix<T> A_new_sub(A.rows - r, A.cols - r);
        for (int i = r; i < A.rows; ++i) {
            for (int j = r; j < A.cols; ++j) {
                A_new_sub.data[i - r][j - r] = A.data[i][j];
            }
        }
        A_new_sub = H * A_new_sub;
        for (int i = r; i < A.rows; ++i) {
            for (int j = r; j < A.cols; ++j) {
                A.data[i][j] = A_new_sub.data[i - r][j - r];
            }
        }
        for (int i = r + 1; i < A.cols; ++i) {
            A.data[r][i] = v.data[2 + r];
        }
        for (int i = r + 1; i < A.rows; ++i) {
            c[i] = c[i] - A.data[i][r] * A.data[i][r];
        }
        t = -MAXFLOAT;
        for (int i = r + 1; i < c.size(); ++i) {
            if (c[i] > t) {
                t = c[i];
            }
        }
    }



}

//helper
template <typename T>
bool has_Converged(Vector<T>& diagonal, const Vector<T>& prevDiagonal, double tolerance) {
    for (size_t i = 0; i < diagonal.size; ++i) {
        if (fabs(diagonal.data.at(i) - prevDiagonal.data.at(i)) > tolerance) {
            return false;
        }
    }
    return true;
}

template <typename T>
Vector<T> back_substitution(Matrix<T> &R, Vector<T> &b) {
    int n = R.rows;
    Vector<T> x(n);

    for (int i = n - 1; i >= 0; --i) {
        T sum = 0;
        for (int j = i + 1; j < n; ++j) {
            sum += R.data[i][j] * x.data[j];
        }
        x.data[i] = (b.data[i] - sum) / R.data[i][i];
    }

    return x;
}

template <typename T>
void Matrix<T>::eigen(Matrix<T> &eigen_vectors, Vector<T> &eigen_values, double tolerance){
    int n = this->cols; 
    Matrix <T> A = *this;

    eigen_values = Vector<T>(n);
    eigen_vectors = Matrix<T>::Identity(n);
    Vector<T> prev_diag(n);

    Matrix<T> Q;
    Matrix<T> R;
    do {
        prev_diag = eigen_values;

        QR_decomp(A, Q, R);
        A = Q.transpose() * A * Q;

        for (int i = 0; i < n; ++i) {
            eigen_values.data.at(i) = A.data.at(i).at(i);
        }
    } while (!has_Converged(eigen_values, prev_diag, tolerance));

}

template <typename T>
Matrix<T> Matrix<T>::singular_value_decomposition(Matrix<T>& U, Matrix<T>& A, Matrix<T>& V) {
    A.print();
    Matrix<T> c = A.transpose() * A;
    V = Matrix<T>::Identity(c.rows);
    for (int i = 0; i < 20; ++i) {
        Matrix<T> Q, R;
        QR_decomp(c, Q, R);
        c = R * Q;
        V = V * Q;
    }
        printf("here \n");
    Matrix<T> AV = A * V;
    Matrix<T> Ut, A_M, Sigma;
    AV.QR_column_pivoting(Ut, A_M, Sigma);
    U = Ut.transpose();
    return Sigma;

}

template <typename T>
Matrix<T> Matrix<T>::moore_penrose() {
    Matrix<T> V, S, U;
    this->singular_value_decomposition(U, S, V);
    Matrix<T> S_plus = S;

    for (int i = 0; i < S_plus.rows; i++) {
        for (int j = 0; j < S_plus.cols; j++) {
            if (S_plus.data.at(i).at(j) != 0) {
                S_plus.data.at(i).at(j) = 1 / S_plus.data.at(i).at(j);
            }
        }
    }
    S_plus = S_plus.transpose();

    Matrix<T> inv;
    inv = V * S_plus * U.transpose();

    return inv;
}

template <typename T>
Matrix<T> Matrix<T>::cholesky() const {
    if (rows != cols) {
        throw std::logic_error("Matrix must be square for Cholesky decomposition");
    }
    Matrix<T> L(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j <= i; ++j) {
            T sum = 0;
            for (int k = 0; k < j; ++k) {
                sum += L.data[i][k] * L.data[j][k];
            }
            if (i == j) {
                L.data[i][j] = sqrt(data[i][i] - sum);
            } else {
                L.data[i][j] = (data[i][j] - sum) / L.data[j][j];
            }
        }
    }

    return L;
}
