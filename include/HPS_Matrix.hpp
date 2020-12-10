#ifndef HPS_MATRIX_HPP_
#define HPS_MATRIX_HPP_

#include <cassert>
#include <initializer_list>
#include <iostream>
#include <random>
#include <type_traits>
#include "HPS_Vector.hpp"

extern "C" {
    void dgemv_(char* TRANS, int* M, int* N, double* ALPHA, double* A, int* LDA, double* X, int* INCX, double* BETA, double* Y, int* INCY);
    void dgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* ALPHA, double* A, int* LDA, double* B, int* LDB, double* BETA, double* C, int* LDC);
    void dgesv_(int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
}

namespace hps {

template <class T>
class Matrix {

public:

    hps::Vector<T> data_{};

    Matrix(std::size_t num_rows, std::size_t num_cols)
        : rows_{num_rows}, cols_{num_cols}, data_(num_rows*num_cols)
            {}

    Matrix(std::size_t num_rows, std::size_t num_cols, T& data)
        : rows_{num_rows}, cols_{num_cols}, data_(num_rows*num_cols, data)
            {}

    Matrix(std::size_t num_rows, std::size_t num_cols, const hps::Vector<T>& data)
        : rows_{num_rows}, cols_{num_cols}, data_{data}
            {}

    Matrix(const Matrix<T>& mat)
        : rows_{mat.rows_}, cols_{mat.cols_}, data_{mat.data_}
            {}

    /**
     * Accesses the element at the given position.
     * @param  row_index Row index
     * @param  col_index Column index
     * @return           Reference to element at indices.
     */
    T& at(std::size_t row_index, std::size_t col_index) noexcept {
        return data_[flattenIndex(row_index, col_index)];
    }

    /**
     * \see at
     */
    const T& at(std::size_t row_index, std::size_t col_index) const noexcept {
        return data_[flattenIndex(row_index, col_index)];
    }

    /**
     * Returns a direct pointer to the underlying Vector's memory array (calls std::vector.data())
     * @return Pointer to memory array
     */
    T* data() noexcept {
        return data_.data();
    }

    /**
     * \see data
     */
    const T* data() const noexcept {
        return data_.data();
    }

    /**
     * Transposes the matrix and returns a new matrix object.
     * @return New transposed matrix object
     */
    Matrix<T> transpose() const {
        hps::Vector<T> data;
        data.reserve(cols_*rows_);

        for (std::size_t j = 0; j < cols_; j++) {
            for (std::size_t i = 0; i < rows_; i++) {
                data.emplace_back(at(i,j));
            }
        }

        return Matrix<T>{cols_, rows_, std::move(data)};
    }

    /**
     * Extracts a block matrix from the matrix.
     * @param  row_index  Starting row index of block to extract
     * @param  col_index  Starting column index of block to extract
     * @param  row_length Number of rows to extract (and number of rows in extracted block matrix)
     * @param  col_length Number of columns to extract (and number of columns in extracted block matrix)
     * @return            New matrix object of extracted block matrix
     */
    Matrix<T> extract(std::size_t row_index, std::size_t col_index, std::size_t row_length, std::size_t col_length) const {
        if (row_index + row_length > rows_) {
            throw std::invalid_argument("[Matrix<T>::extract] Row size exceeds matrix size");
        }
        if (col_index + col_length > cols_) {
            throw std::invalid_argument("[Matrix<T>::extract] Column size exceeds matrix size");
        }

        hps::Vector<T> out_data;
        out_data.reserve(row_length * col_length);
        for (std::size_t i = row_index; i < row_index + row_length; i++) {
            for (std::size_t j = col_index; j < col_index + col_length; j++) {
                out_data.emplace_back(at(i,j));
            }
        }

        return {row_length, col_length, std::move(out_data)};
    }

    /**
     * Opposite of extract. Puts a matrix into `this` matrix at the specified row and column indices
     * @param row_index Row index of location to insert matrix
     * @param col_index Column index of location to insert matrix
     * @param mat       Matrix to insert
     */
    void intract(std::size_t row_index, std::size_t col_index, const Matrix<T>& mat) {
        if (row_index + mat.rows_ > rows_) {
            throw std::invalid_argument("[Matrix<T>::intract] Row size exceeds matrix size");
        }
        if (col_index + mat.cols_ > cols_) {
            throw std::invalid_argument("[Matrix<T>::intract] Column size exceeds matrix size");
        }

        for (std::size_t i = row_index; i < row_index + mat.rows_; i++) {
            for (std::size_t j = col_index; j < col_index + mat.cols_; j++) {
                at(i,j) = mat.at(i - row_index, j - col_index);
            }
        }
    }

    /**
     * Randomizes the data of the matrix.
     * @param lower_bound Lower bound of random data
     * @param upper_bound Upper bound of random data
     */
    void randomize(T lower_bound, T upper_bound) {
        std::mt19937_64 engine{std::random_device{}()};
        std::uniform_real_distribution<T> dist(lower_bound, upper_bound);
        for(std::size_t i = 0; i < rows_; i++) {
            for(std::size_t j = 0; j < cols_; j++) {
                data_[flattenIndex(i, j)] = dist(engine);
            }
        }
    }

    /**
     * Negates all data.
     */
    void negate() {
        for (std::size_t i = 0; i < data_.size(); i++) {
            data_[i] = -data_[i];
        }
    }

    /**
     * Sets the matrix to be the Identity matrix.
     */
    void identify() {
        for (std::size_t i = 0; i < rows_; i++) {
            for (std::size_t j = 0; j < cols_; j++) {
                if (i == j) {
                    this->at(i,j) = 1;
                }
                else {
                    this->at(i,j) = 0;
                }
            }
        }
    }

    /**
     * Returns the largest value in the matrix.
     * @return Largest value in matrix
     */
    T max() {
        T m = this->data_[0];
        for (std::size_t i = 0; i < this->size(); i++) {
            if (this->data_[i] > m) {
                m = this->data_[i];
            }
        }
        return m;
    }

    /**
     * Returns the number of rows in the matrix.
     * @return Number of rows
     */
    std::size_t rows() const noexcept {
        return rows_;
    }

    /**
     * Returns the number of columns in the matrix.
     * @return Number of columns
     */
    std::size_t cols() const noexcept {
        return cols_;
    }

    /**
     * Returns the number of elements in the matrix.
     * @return Number of elements in the matrix
     */
    std::size_t size() const noexcept {
        return data_.size();
    }

    /**
     * Copy assignment operator.
     * @param other Matrix to assign from
     * @return      Self after assigning from other
     */
    Matrix<T>& operator=(const Matrix<T>& other) {
        rows_ = other.rows_;
        cols_ = other.cols_;
        data_ = other.data_;
        return *this;
    }

    /**
     * Move assignment operator. Note that other has a size of zero after being move from.
     * @param other Matrix to assign from
     * @return      Self after assigning from other
     */
    Matrix<T>& operator=(const Matrix<T>&& other) noexcept {
        rows_ = other.rows_;
        cols_ = other.cols_;
        data_ = std::move(other.data_);
        // other.rows_ = 0;
        // other.cols_ = 0;
        return *this;
    }

    /**
     * Tests equality of `this` and other matrix
     * @param  other Other matrix to compare to `this`
     * @retval true  IFF matrices are the same size and have the same data
     * @retval false Otherwise
     */
    bool operator==(const Matrix<T>& other) const {
        if (rows_ != other.rows_) return false;
        if (cols_ != other.cols_) return false;
        if (data_.size() != other.data_.size()) return false;

        for (std::size_t i = 0; i < data_.size(); i++) {
            if (data_[i] != other.data_[i]) return false;
        }

        return true;
    }

    /**
     * Inverse of == operator. \see operator==
     */
    bool operator!=(const Matrix<T>& other) const {
        return !(*this == other);
    }

    /**
     * Adds supplied matrix to this. Matrices must have the same dimensions.
     * @param  rhs Matrix to add to `this`
     * @return     `this` after adding `rhs` to it
     */
    Matrix<T>& operator+=(const Matrix<T>& rhs) {
        if (rows_ != rhs.rows_ || cols_ != rhs.cols_ || data_.size() != rhs.data_.size()) {
            throw std::invalid_argument("[Matrix<T>::operator+=] Matrices must be the same size");
        }

        for (std::size_t i = 0; i < data_.size(); i++) {
            data_[i] += rhs.data_[i];
        }
        return *this;
    }

    /**
     * Subtracts supplied matrix to `this`. Matrices must have the same dimensions.
     * @param  rhs Matrix to subtract from `this`
     * @return     `this` after subtracting `rhs` from it
     */
    Matrix<T>& operator-=(const Matrix<T>& rhs) {
        if (rows_ != rhs.rows_ || cols_ != rhs.cols_ || data_.size() != rhs.data_.size()) {
            throw std::invalid_argument("[Matrix<T>::operator-=] Matrices must be the same size");
        }

        for (std::size_t i = 0; i < data_.size(); i++) {
            data_[i] -= rhs.data_[i];
        }
        return *this;
    }

    /**
     * Basic printing to streams
     * @param  os  Stream to print to
     * @param  mat Matrix to print
     * @return     Stream after printing matrix
     */
    friend std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat) {
        os << '{';

        for (std::size_t i = 0; i < mat.rows_; i++){
            for (std::size_t j = 0; j < mat.cols_; j++){
                os << mat.at(i, j);
                os << ',';
                // if(j != mat.cols_ - 1){
                //     os << ',';
                // }
            }
            if(i != mat.rows_-1){
                os << '\n';
            }
        }

        os << '}';
        return os;
    }

private:

    std::size_t rows_{};
    std::size_t cols_{};

    /**
     * Flattens the two dimensional index to a one dimensional index. Asserts that row_index and col_index are within bounds and asserts that the flattened index is also within bounds.
     * @param  row_index Row index
     * @param  col_index Column index
     * @return           Flattened index for 1D data access
     */
    std::size_t flattenIndex(std::size_t row_index, std::size_t col_index) const noexcept {
        assert(row_index < rows_);
        assert(col_index < cols_);

        const std::size_t index = row_index*cols_ + col_index;
        assert(index < rows_*cols_);

        return index;
    }

};

// @@TODO: Add +, -, and * operators
/**
 * Adds two matrices.
 * \see Matrix<T>::operator+=
 */
template<class T>
Matrix<T> operator+(Matrix<T> lhs, const Matrix<T>& rhs) {
    lhs += rhs;
    return lhs;
}

/**
 * Subtracts two matrices.
 * \see Matrix<T>::operator-=
 */
template<class T>
Matrix<T> operator-(Matrix<T> lhs, const Matrix<T>& rhs) {
    lhs -= rhs;
    return lhs;
}

/**
 * Multiplies a matrix `A` and a vector `x`. Wrapper for the BLAS dgemv.
 * @param   A  Matrix to multiply
 * @param   x  Vector to multiply
 * @return     New vector instance with result
 */
template<class T>
hps::Vector<T> operator*(Matrix<T>& A, hps::Vector<T>& x) {

    // Check inputs
    if (A.cols() != x.size()) {
        throw std::invalid_argument("[hps::Vector<T> operator*] Invalid matrix and vector dimensions. Matrix `A` must have the same number of columns as the size of Vector `x`");
    }

    // Create output vector
    hps::Vector<T> b(A.rows());

    // Setup call
    // void dgemv_(int* TRANS, int* M, int* N, double* ALPHA, double* A, int* LDA, double* X, int* INCX, double* BETA, double* Y, int* INCY);
    char TRANS_ = 'C';
    int M_ = A.cols();
    int N_ = A.rows();
    double ALPHA_ = 1.0;
    double* A_ = A.data();
    int LDA_ = A.cols();
    double* X_ = x.data();
    int INCX_ = 1;
    double BETA_ = 0.0;
    double* Y_ = b.data();
    int INCY_ = 1;
    dgemv_(&TRANS_, &M_, &N_, &ALPHA_, A_, &LDA_, X_, &INCX_, &BETA_, Y_, &INCY_);

    return b;
}

/**
 * Multiplies a matrix and a matrix. Wrapper for LAPACK dgemm.
 * @param  A  LHS matrix to multiply
 * @param  B  RHS matrix to multiply
 * @return    New matrix instance with result
 */
template<class T>
Matrix<T> operator*(Matrix<T>& A, Matrix<T>& B) {

    // Check inputs
    if (A.cols() != B.rows()) {
        throw std::invalid_argument("[hps::Matrix<T> operator*] Invalid matrix dimensions. Matrix `A` must have the same number of columns as the rows of Matrix `B`");
    }

    // Create output vector
    Matrix<T> C(A.rows(), B.cols());

    // Setup call
    // void dgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* ALPHA, double* A, int* LDA, double* B, int* LDB, double* BETA, dobule* C, int* LDC);
    char TRANSA_ = 'C';
    char TRANSB_ = 'C';
    int M_ = C.cols();
    int N_ = C.cols();
    int K_ = B.rows();
    double ALPHA_ = 1.0;
    double* A_ = A.data();
    int LDA_ = K_;
    double* B_ = B.data();
    int LDB_ = N_;
    double BETA_ = 0.0;
    double* C_ = C.data();
    int LDC_ = M_;
    dgemm_(&TRANSA_, &TRANSB_, &M_, &N_, &K_, &ALPHA_, A_, &LDA_, B_, &LDB_, &BETA_, C_, &LDC_);

    return C.transpose();
}

template<class T>
hps::Vector<T> solve(Matrix<T>& A, hps::Vector<T>& b) {

    // Check inputs
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("[hps::Matrix<T> solve] Matrix must be square (rows != cols).");
    }

    // Setup output vector and utility variables
    hps::Vector<T> x(b);
    Matrix<T> AT = A.transpose();
    Vector<int> p(b.size());

    // Setup call
    int N_ = AT.rows();
    int NRHS_ = 1;
    double* A_ = AT.data();
    int LDA_ = AT.rows();
    int* IPIV_ = p.data();
    double* B_ = x.data();
    int LDB_ = b.size();
    int INFO_;
    dgesv_(&N_, &NRHS_, A_, &LDA_, IPIV_, B_, &LDB_, &INFO_);

    // Check output
    if (INFO_ != 0) {
        std::cerr << "[hps::Vector<T> solve] Fortran call to `dgesv_` returned non-zero flag of: " << INFO_ << std::endl;
    }

    return x;
}

template<class T>
Matrix<T> solve(Matrix<T>& A, Matrix<T>& B) {

    // Check inputs
    if (A.rows() != A.cols()) {
        throw std::invalid_argument("[hps::Matrix<T> solve] Matrix must be square (rows != cols).");
    }

    // Setup output matrix and utility variables
    hps::Matrix<T> X = B.transpose();
    hps::Matrix<T> AT = A.transpose();
    hps::Vector<int> p(A.rows());

    // Setup call
    int N_ = AT.rows();
    int NRHS_ = B.cols();
    double* A_ = AT.data();
    int LDA_ = AT.rows();
    int* IPIV_ = p.data();
    double* B_ = X.data();
    int LDB_ = B.rows();
    int INFO_;
    dgesv_(&N_, &NRHS_, A_, &LDA_, IPIV_, B_, &LDB_, &INFO_);

    // Check output
    if (INFO_ != 0) {
        std::cerr << "[hps::Matrix<T> solve] Fortran call to `dgesv_` returned non-zero flag of: " << INFO_ << std::endl;
    }

    return X.transpose();
}

}

#endif // HPS_MATRIX_HPP_
