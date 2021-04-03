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
add_subdirectory(external)

# C++17 Standard
set(CMAKE_CXX_STANDARD 17)
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
include_directories("${PROJECT_SOURCE_DIR}/external")

# Library Directories
link_directories(${MY_LIB_DIR})

# Create Source Files
file(GLOB HPS_SOURCES "${PROJECT_SOURCE_DIR}/src/HPS_*.cpp")

# Create App Executables
#    Convergence Analysis
add_executable(convergence ${HPS_SOURCES} ${APP_DIR}/convergence.cpp)
target_link_libraries(convergence gfortran fish90 blas lapack lmmio)

#    4-to-1 Merge (DtN check)
add_executable(4to1_merge ${HPS_SOURCES} ${APP_DIR}/4to1_merge.cpp)
target_link_libraries(4to1_merge gfortran fishpack blas lapack lmmio)

#    HPS 4-to-1 build and solve
add_executable(HPS ${HPS_SOURCES} ${APP_DIR}/HPS.cpp)
target_link_libraries(HPS gfortran fishpack blas lapack lmmio)

#    Main Driver (Testing and Experimenting right now)
add_executable(main ${HPS_SOURCES} ${APP_DIR}/main.cpp)
target_link_libraries(main gfortran fish90 blas lapack lmmio)
