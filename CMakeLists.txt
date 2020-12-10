# CMake Version
cmake_minimum_required(VERSION 3.18)

# Project Info
project(HPS VERSION 0.1 LANGUAGES C CXX Fortran)

# Project Variables
set(APP_DIR src/apps)
set(MY_LIB_DIR "/home/dchipman1/.local/lib")

# External Packages
include(FindBLAS)
include(FindLAPACK)

# C++11 Standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Configure GTest
configure_file(CMakeLists.txt.in GTestFiles/googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/GTestFiles/googletest-download )
execute_process(COMMAND ${CMAKE_COMMAND} --build . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/GTestFiles/googletest-download )
add_subdirectory(${CMAKE_BINARY_DIR}/GTestFiles/googletest-src ${CMAKE_BINARY_DIR}/GTestFiles/googletest-build)

# Enable Testing
enable_testing()
add_subdirectory(test)

# Build File Structure
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Include Directories
include_directories("${PROJECT_SOURCE_DIR}/include")

# Library Directories
link_directories(${MY_LIB_DIR})

# Create Source Files
file(GLOB HPS_SOURCES "${PROJECT_SOURCE_DIR}/src/HPS_*.cpp")

# Create App Executables
#    Convergence Analysis
add_executable(convergence ${HPS_SOURCES} ${APP_DIR}/convergence.cpp)
target_link_libraries(convergence gfortran fish90 blas lapack)

#    Main Driver (Testing and Experimenting right now)
add_executable(main ${HPS_SOURCES} ${APP_DIR}/main.cpp)
target_link_libraries(main gfortran fish90 blas lapack)
