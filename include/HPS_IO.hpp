#ifndef HPS_IO_HPP_
#define HPS_IO_HPP_

#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdio>

extern "C"{
#include "mmio/include/mmio.h"
}

#include "HPS_Matrix.hpp"
#include "HPS_Vector.hpp"
#include "HPS_Grid.hpp"

namespace hps {

enum class DataType {integer, real};

//0 if file is valid
//1 if file cannot be opened
//2 if file doesn't contain valid banner
//3 if file attributes are not supported.
int validFile(const std::string& fileName);

DataType readDataType(const std::string& fileName);

Matrix<int> readMatrixInt(std::FILE* file);
Matrix<double> readMatrixDouble(std::FILE* file);

bool writeMatrixInt(const std::string& fileName, const Matrix<int>& mat);
bool writeMatrixDouble(const std::string& fileName, const Matrix<double>& mat);

Vector<int> readVectorInt(std::FILE* file);
Vector<double> readVectorDouble(std::FILE* file);

bool writeVectorInt(const std::string& fileName, const Vector<int>& vec);
bool writeVectorDouble(const std::string& fileName, const Vector<double>& vec);

bool writeVTK(CellGrid<double, 2> grid, Matrix<double> matrix, std::string& fileName, std::string& grid_name, std::string& data_name);

} // END NAMESPACE hps

#endif // HPS_IO_HPP_
