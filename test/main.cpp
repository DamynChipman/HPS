#include "gtest/gtest.h"

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    // ::testing::internal::CaptureStdout();
    return RUN_ALL_TESTS();
}
