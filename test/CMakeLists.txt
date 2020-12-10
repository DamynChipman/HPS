include_directories("${PROJECT_SOURCE_DIR}/include")

file(GLOB HPS_SOURCES "${PROJECT_SOURCE_DIR}/src/HPS_*.cpp")
file(GLOB TESTS "${PROJECT_SOURCE_DIR}/test/Test*.cpp")

# External Packages
include(FindBLAS)
include(FindLAPACK)

set(MY_LIB_DIR "/home/dchipman1/.local/lib")
link_directories(${MY_LIB_DIR})

include(GoogleTest)

foreach(file ${TESTS})
  set(name)
  get_filename_component(name ${file} NAME_WE)
  add_executable("${name}_tests"
    ${HPS_SOURCES}
    ${file}
    "${PROJECT_SOURCE_DIR}/test/main.cpp")
  target_link_libraries("${name}_tests" gtest_main fish90 gfortran blas lapack)
  gtest_discover_tests("${name}_tests")
  # add_test(NAME ${name} COMMAND "${name}_tests")
endforeach()
