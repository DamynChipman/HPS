#ifndef HPS_VECTOR_HPP_
#define HPS_VECTOR_HPP_

#include <vector>
#include <iostream>
#include <cmath>

namespace hps {

template <class T>
class Vector : public std::vector<T> {

public:

    using std::vector<T>::vector;

    Vector<T> extract(std::size_t start_index, std::size_t length) const {
        if (start_index + length > this->std::vector<T>::size()) {
            throw std::invalid_argument("[Vector<T>::extract] Size of vector to extract exceeds original vector size");
        }

        Vector<T> out_vector;
        out_vector.reserve(length);
        for (std::size_t i = start_index; i < start_index + length; i++) {
            out_vector.emplace_back(this->std::vector<T>::operator[](i));
        }

        return out_vector;
    }

    void intract(std::size_t start_index, const Vector<T>& vec) {
        if (start_index + vec.size() > this->std::vector<T>::size()) {
            throw std::invalid_argument("[Vector<T>::intract] Vector to intract size is larger than size of host vector");
        }

        for (std::size_t i = start_index; i < start_index + vec.size(); i++) {
            this->std::vector<T>::operator[](i) = vec[i - start_index];
        }
    }

    T max() const {
        T to_return;
        if (this->std::vector<T>::size() == 1) {
            to_return = this->std::vector<T>::operator[](0);
        }
        else {
            to_return = this->std::vector<T>::operator[](0);
            for (std::size_t i = 0; i < std::vector<T>::size(); i++) {
                if (std::vector<T>::operator[](i) > to_return) {
                    to_return = std::vector<T>::operator[](i);
                }
            }
        }
        return to_return;
    }

    T min() const {
        T to_return = this->std::vector<T>::operator[](0);;
        if (this->std::vector<T>::size() == 1) {
            return to_return;
        }
        else {
            for (std::size_t i = 0; i < std::vector<T>::size(); i++) {
                if (std::vector<T>::operator[](i) < to_return) {
                    to_return = std::vector<T>::operator[](i);
                }
            }
        }
        return to_return;
    }

    Vector<T> abs() {
        Vector<T> to_return(*this);
        for (std::size_t i = 0; i < std::vector<T>::size(); i++) {
            to_return[i] = fabs(std::vector<T>::operator[](i));
        }
        return to_return;
    }

    double norm() {
        T to_return = 0;
        for (std::size_t i = 0; i < std::vector<T>::size(); i++) {
            to_return += pow(fabs(std::vector<T>::operator[](i)), 2);
        }
        return sqrt(to_return);
    }

    double gridNorm(double delta_x, double delta_y, int order) {
        if (order <= 0) {
            throw std::invalid_argument("[Vector<T>::gridNorm] Invalid grid norm `order`. Options are a positive, non-zero integer or `inifity`");
        }

        double to_return = 0;
        for (std::size_t i = 0; i < std::vector<T>::size(); i++) {
            to_return += pow(fabs(std::vector<T>::operator[](i)), order);
        }
        to_return = pow(delta_x * delta_y * to_return, (double) 1/order);
        return to_return;
    }

    double gridNorm(double delta_, double delta_y, std::string order) {
        if (order != "infinity") {
            throw std::invalid_argument("[Vector<T>::gridNorm] Invalid grid norm `order`. Options are a positive, non-zero integer or `inifity`");
        }

        return this->abs().max();
    }

    Vector<T>& operator+=(const Vector<T>& rhs) {
        if (this->std::vector<T>::size() != rhs.size()) {
            throw std::invalid_argument("[Vector<T>::operator+=] Vectors must be the same size");
        }

        for (std::size_t i = 0; i < this->std::vector<T>::size(); i++) {
            this->std::vector<T>::operator[](i) += rhs[i];
        }
        return *this;
    }

    Vector<T>& operator-=(const Vector<T>& rhs) {
        if (this->std::vector<T>::size() != rhs.size()) {
            throw std::invalid_argument("[Vector<T>::operator+=] Vectors must be the same size");
        }

        for (std::size_t i = 0; i < this->std::vector<T>::size(); i++) {
            this->std::vector<T>::operator[](i) -= rhs[i];
        }
        return *this;
    }

    Vector<T>& operator*=(const double rhs) {
        for (std::size_t i = 0; i < this->std::vector<T>::size(); i++) {
            this->std::vector<T>::operator[](i) *= rhs;
        }
        return *this;
    }

    Vector<T>& operator*=(const int rhs) {
        for (std::size_t i = 0; i < this->std::vector<T>::size(); i++) {
            this->std::vector<T>::operator[](i) *= rhs;
        }
        return *this;
    }

    friend std::ostream& operator<<(std:: ostream& os, const Vector<T>& vec) {
        os << "{";
        for (std::size_t i = 0; i < vec.size(); i++) {
            os << vec[i];
            if (i != vec.size() - 1) {
                os << ", ";
            }
        }
        os << "}";
        return os;
    }

};

template<class T>
Vector<T> operator+(Vector<T> lhs, const Vector<T>& rhs) {
    lhs += rhs;
    return lhs;
}

template<class T>
Vector<T> operator-(Vector<T> lhs, const Vector<T>& rhs) {
    lhs -= rhs;
    return lhs;
}

template<class T>
Vector<T> operator*(Vector<T> lhs, const double rhs) {
    lhs *= rhs;
    return lhs;
}

template<class T>
Vector<T> operator*(Vector<T> lhs, const int rhs) {
    lhs *= rhs;
    return lhs;
}

}

#endif // HPS_VECTOR_HPP_
