#include "gtest/gtest.h"
#include "HPS_Vector.hpp"

#include <iostream>
#include <string>

/*///////////////////////////////////////////////
 *
 *  Unit Tests -- Vector Class
 *
 *///////////////////////////////////////////////

TEST(Vector, size) {
    std::size_t n = 5;
    hps::Vector<int> vec(n);
    ASSERT_EQ(vec.size(), n);
}

TEST(Vector, const_data) {
    std::size_t n = 5;
    int data = 42;
    hps::Vector<int> vec(n, 42);
    for (std::size_t i = 0; i < n; i++) {
        ASSERT_EQ(vec[i], data);
    }
}

TEST(Vector, variable_data) {
    std::size_t n = 5;
    hps::Vector<int> vec({0,1,2,3,4});
    for (int i = 0; i < n; i++) {
        ASSERT_EQ(vec[i], i);
    }
}

TEST(Vector, extract) {
    std::size_t n = 6;
    hps::Vector<int> vec({0,1,2,3,4,5});
    hps::Vector<int> extracted = vec.extract(2, 2);
    hps::Vector<int> expected({2,3});
    
    for (int i = 0; i < 2; i++) {
        ASSERT_EQ(extracted[i], expected[i]);
    }
}

TEST(Vector, intract) {
    std::size_t n = 5;
    hps::Vector<int> vec({0,1,2,3,4});
    hps::Vector<int> to_intract({10,20,30});
    hps::Vector<int> expected({0,10,20,30,4});
    
    vec.intract(1, to_intract);
    
    for (int i = 0; i < 5; i++) {
        ASSERT_EQ(vec[i], expected[i]);
    }
}

TEST(Vector, vec_plus_vec) {
    hps::Vector<int> lhs({0,1,2,3});
    hps::Vector<int> rhs({0,10,20,30});
    hps::Vector<int> expected({0,11,22,33});
    hps::Vector<int> test = lhs + rhs;
    for (int i = 0; i < 4; i++) {
        ASSERT_EQ(test[i], expected[i]);
    }
}

TEST(Vector, vec_minus_vec) {
    hps::Vector<int> lhs({0,1,2,3});
    hps::Vector<int> rhs({0,10,20,30});
    hps::Vector<int> expected({0,-9,-18,-27});
    hps::Vector<int> test = lhs - rhs;
    for (int i = 0; i < 4; i++) {
        ASSERT_EQ(test[i], expected[i]);
    }
}

TEST(Vector, vec_times_intvalue) {
    hps::Vector<int> lhs({0,1,2,3});
    int rhs = 2;
    hps::Vector<int> expected({0,2,4,6});
    lhs *= rhs;
    for (int i = 0; i < 4; i++) {
        ASSERT_EQ(lhs[i], expected[i]);
    }
}

TEST(Vector, vec_times_doublevalue) {
    hps::Vector<double> lhs({0,1,2,3});
    double rhs = 2;
    hps::Vector<double> expected({0,2,4,6});
    lhs *= rhs;
    for (int i = 0; i < 4; i++) {
        ASSERT_EQ(lhs[i], expected[i]);
    }
}

TEST(Vector, max) {
    hps::Vector<double> vec({-0.1, 2.4, 0.0, -9.8, 42.0});
    double max_true = 42.0;
    double max_test = vec.max();
    ASSERT_EQ(max_test, max_true);
}