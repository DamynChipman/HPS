#include "gtest/gtest.h"
#include "HPS_Matrix.hpp"
#include "HPS_Vector.hpp"
#include "HPS_Patch.hpp"

#include <iostream>
#include <string>

/*///////////////////////////////////////////////
 *
 *  Unit Tests -- Matrix Class
 *
 *///////////////////////////////////////////////

TEST(Matrix, sizes) {
    std::size_t r = 2;
    std::size_t c = 3;
    hps::Matrix<int> mat(r, c);

    ASSERT_EQ(mat.rows(), r);
    ASSERT_EQ(mat.cols(), c);
}

TEST(Matrix, const_data) {
    std::size_t r = 3;
    std::size_t c = 4;
    int data = 42;
    hps::Matrix<int> mat(r, c, data);

    for (std::size_t i = 0; i < r; i++) {
        for (std::size_t j = 0; j < c; j++) {
            ASSERT_EQ(mat.at(i,j), data);
        }
    }
}

TEST(Matrix, variable_data) {
    std::size_t r = 3;
    std::size_t c = 4;
    hps::Matrix<int> mat(r, c,
        {1, 2,  3,  4,
         5, 6,  7,  8,
         9, 10, 11, 12});

    int data = 1;
    for (std::size_t i = 0; i < r; i++) {
        for (std::size_t j = 0; j < c; j++) {
            ASSERT_EQ(mat.at(i,j), data);
            data++;
        }
    }
}

TEST(Matrix, transpose_size) {
    std::size_t r = 3;
    std::size_t c = 4;
    hps::Matrix<int> mat(r, c,
        {1, 2,  3,  4,
         5, 6,  7,  8,
         9, 10, 11, 12});

    ASSERT_EQ(mat.transpose().rows(), c);
    ASSERT_EQ(mat.transpose().cols(), r);
}

TEST(Matrix, transpose_data) {
    std::size_t r = 3;
    std::size_t c = 4;
    hps::Matrix<int> mat(r, c,
        {1, 2,  3,  4,
         5, 6,  7,  8,
         9, 10, 11, 12});

    hps::Matrix<int> mat_T(c, r,
        {1,  5,  9,
         2,  6,  10,
         3,  7,  11,
         4,  8,  12});

    for (std::size_t i = 0; i < c; i++) {
        for (std::size_t j = 0; j < r; j++) {
            ASSERT_EQ(mat.transpose().at(i,j), mat_T.at(i,j));
        }
    }
}

TEST(Matrix, is_equal_operator) {
    std::size_t r = 3;
    std::size_t c = 4;
    hps::Matrix<int> mat1(r, c,
        {1, 2,  3,  4,
         5, 6,  7,  8,
         9, 10, 11, 12});
    hps::Matrix<int> mat2(r, c,
        {1, 2,  3,  4,
         5, 6,  7,  8,
         9, 10, 11, 12});

    ASSERT_EQ(mat1, mat2);
}

TEST(Matrix, identify) {
    int r = 3;
    int c = 3;
    hps::Matrix<int> mat1(r, c);
    hps::Matrix<int> mat2(r, c,
        {1, 0, 0,
         0, 1, 0,
         0, 0, 1});

    mat1.identify();
    ASSERT_EQ(mat1, mat2);
}

TEST(Matrix, negate) {
    int r = 3;
    int c = 4;
    hps::Matrix<int> mat1(r, c,
        {1, 2,  3,  4,
         5, 6,  7,  8,
         9, 10, 11, 12});
    hps::Matrix<int> mat2(r, c,
        {-1, -2,  -3,  -4,
         -5, -6,  -7,  -8,
         -9, -10, -11, -12});

    mat1.negate();
    ASSERT_EQ(mat1, mat2);
}

TEST(Matrix, extract_block) {
    std::size_t r = 3;
    std::size_t c = 4;
    hps::Matrix<int> mat(r, c,
        {1, 2,  3,  4,
         5, 6,  7,  8,
         9, 10, 11, 12});

    hps::Matrix<int> mat_1_test = mat.extract(0, 0, 2, 2);
    hps::Matrix<int> mat_2_test = mat.extract(0, 2, 2, 2);
    hps::Matrix<int> mat_3_test = mat.extract(2, 0, 1, 4);

    hps::Matrix<int> mat_1_true(2, 2,
        {1, 2,
         5, 6});
    hps::Matrix<int> mat_2_true(2, 2,
        {3, 4,
         7, 8});
    hps::Matrix<int> mat_3_true(1, 4,
        {9, 10, 11, 12});

    ASSERT_EQ(mat_1_test, mat_1_true);
    ASSERT_EQ(mat_2_test, mat_2_true);
    ASSERT_EQ(mat_3_test, mat_3_test);
}

TEST(Matrix, intract_block) {
    std::size_t r = 3;
    std::size_t c = 4;
    hps::Matrix<int> mat_test(r, c);

    hps::Matrix<int> mat_1(2, 2,
        {1, 2,
         5, 6});
    hps::Matrix<int> mat_2(2, 2,
        {3, 4,
         7, 8});
    hps::Matrix<int> mat_3(1, 4,
        {9, 10, 11, 12});

    mat_test.intract(0, 0, mat_1);
    mat_test.intract(0, 2, mat_2);
    mat_test.intract(2, 0, mat_3);

    hps::Matrix<int> mat_true(r, c,
        {1, 2,  3,  4,
         5, 6,  7,  8,
         9, 10, 11, 12});

    ASSERT_EQ(mat_test, mat_true);
}

TEST(Matrix, add) {
    int r = 3;
    int c = 4;
    hps::Matrix<int> mat(r, c,
        {1, 2,  3,  4,
         5, 6,  7,  8,
         9, 10, 11, 12});

    hps::Matrix<int> mat_true(r, c,
        {2,  4,  6,  8,
         10, 12, 14, 16,
         18, 20, 22, 24});

    ASSERT_EQ(mat + mat, mat_true);
}

TEST(Matrix, subtract) {
    int r = 3;
    int c = 4;
    hps::Matrix<int> mat(r, c,
        {1, 2,  3,  4,
         5, 6,  7,  8,
         9, 10, 11, 12});

    hps::Matrix<int> mat_true(r, c,
        {0,0,0,0,0,0,0,0,0,0,0,0});

    ASSERT_EQ(mat - mat, mat_true);
}

TEST(Matrix, max) {
    int r = 3;
    int c = 4;
    hps::Matrix<int> mat(r, c,
        {1, 2,  3,  4,
         5, 6,  50,  8,
         9, 10, 11, 12});

    ASSERT_EQ(mat.max(), 50);
}

TEST(Matrix, matvec_multiply) {

    int r = 3;
    int c = 4;
    hps::Matrix<double> A(r, c,
        {1, 2,  3,  4,
         5, 6,  7,  8,
         9, 10, 11, 12});

    hps::Vector<double> x({1, 1, 1, 1});
    hps::Vector<double> b_true({10, 26, 42});
    hps::Vector<double> b_test = A*x;

    for (int i = 0; i < 3; i++) {
        ASSERT_FLOAT_EQ(b_test[i], b_true[i]);
    }
}

TEST(Matrix, matvec_multiply2) {
    int r = 3;
    int c = 3;
    hps::Matrix<double> A(r, c);
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            A.at(i,j) = i + j + 0.5;
        }
    }

    hps::Vector<double> x(r);
    for (int i = 0; i < r; i++) {
        x[i] = (double) 2*i + 0.5;
    }

    hps::Vector<double> b_test = A * x;
    hps::Vector<double> b_true(r);
    b_true[0] = 15.25;
    b_true[1] = 22.75;
    b_true[2] = 30.25;

    for (int i = 0; i < r; i++) {
        ASSERT_FLOAT_EQ(b_test[i], b_true[i]);
    }

}

TEST(Matrix, matmat_multiply) {

    int r = 3;
    int c = 4;
    hps::Matrix<double> A(r, c,
        {1, 2,  3,  4,
         5, 6,  7,  8,
         9, 10, 11, 12});

    hps::Matrix<double> B(c, r,
        {1,  2,  3,
         4,  5,  6,
         7,  8,  9,
         10, 11, 12});

    hps::Matrix<double> C_true(r, r,
        {70,  80,  90,
         158, 184, 210,
         246, 288, 330});

    hps::Matrix<double> C_test = A*B;

    for (int i = 0; i < C_test.rows(); i++) {
        for (int j = 0; j < C_test.cols(); j++) {
            ASSERT_FLOAT_EQ(C_test.at(i,j), C_true.at(i,j));
        }
    }
}

TEST(Matrix, matvec_solve) {
    int N = 3;
    hps::Matrix<double> A(N, N,
        {1, 2, 3,
         4, 5, 6,
         7, 8, 10});
    hps::Vector<double> b({1, 1, 1});

    hps::Vector<double> x_test = hps::solve(A, b);
    hps::Vector<double> x_true({-1, 1, 0});

    for (int i = 0; i < N; i++) {
        EXPECT_NEAR(x_test[i], x_true[i], 1e-15);
    }
}

TEST(Matrix, matmat_solve) {
    int N = 5;
    int NRHS = 3;
    hps::Matrix<double> A(N, N,
        {
            6.80, -6.05, -0.45,  8.32, -9.67,
           -2.11, -3.30,  2.58,  2.71, -5.14,
            5.66, 5.36, -2.70,  4.35, -7.26,
            5.97, -4.44,  0.27, -7.17, 6.08,
            8.23, 1.08,  9.04,  2.14, -6.87
        }
    );
    hps::Matrix<double> B(N, NRHS,
        {
            4.02, -1.56, 9.81,
            6.19,  4.00, -4.09,
           -8.22, -8.67, -4.57,
           -7.57,  1.75, -8.61,
           -3.03,  2.86, 8.99
       }
    );
    hps::Matrix<double> X_true(N, NRHS,
        {
            -0.80071403, -0.38962139,  0.95546491,
            -0.69524338, -0.55442713,  0.22065963,
            0.59391499,   0.84222739,  1.90063673,
            1.32172561,  -0.10380185,  5.35766149,
            0.5657562 ,   0.10571095,  4.04060266
       }
    );

    hps::Matrix<double> X_test = hps::solve(A, B);
    for(int i = 0; i < N; i++) {
        for (int j = 0; j < NRHS; j++) {
            EXPECT_NEAR(X_test.at(i,j), X_true.at(i,j), 1e-8);
        }
    }
}
