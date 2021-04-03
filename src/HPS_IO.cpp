#include "HPS_IO.hpp"

namespace hps {

template <class T>
Matrix<T> readMatrix(std::FILE* file, const std::string& formatCode){
    int row = 0;
    int col = 0;

    if(mm_read_mtx_array_size(file, &row, &col)){
        throw std::runtime_error("Error reading array size");
    }

    Matrix<T> mat({static_cast<std::size_t>(row), static_cast<std::size_t>(col)});

    const std::string format = formatCode + "\n";

    //Read as column major
    for(std::size_t x = 0; x < static_cast<std::size_t>(col); ++x){
        for(std::size_t y = 0; y < static_cast<std::size_t>(row); ++y){
            T temp;
            fscanf(file, format.c_str(), &temp);
            mat.at(y, x) = temp;
        }
    }

    return mat;
}

template <typename T>
Vector<T> readVector(std::FILE* file, const std::string& formatCode){
    int vlen = 0;
    int one = 0;   //ignore
    if(mm_read_mtx_array_size(file, &one, &vlen)){
        throw std::runtime_error("Error reading array size");
    }

    const std::size_t len = static_cast<std::size_t>(vlen);

    if(one != 1){
        throw std::invalid_argument("File does not contain a vector! (number of rows is not 1)");
    }

    const std::string format = formatCode + "\n";

    Vector<T> data;
    data.reserve(len);
    for(std::size_t i = 0; i < len; ++i){
        T temp;
        fscanf(file, format.c_str(), &temp);
        data.emplace_back(temp);
    }

    return data;
}


void writeHeader(std::FILE* file, DataType type){
    MM_typecode code{};

    mm_initialize_typecode(&code);
    mm_set_matrix(&code);
    mm_set_array(&code);

    if(type == DataType::integer){
        mm_set_integer(&code);
    }else{
        mm_set_real(&code);
    }

    mm_write_banner(file, code);
}

template <typename T>
void writeMatrix(std::FILE* file, const Matrix<T>& mat, const std::string& formatCode){
    mm_write_mtx_array_size(file, static_cast<int>(mat.rows()), static_cast<int>(mat.cols()));

    const std::string format = formatCode + "\n";

    for(std::size_t x = 0; x < mat.cols(); ++x){
        for(std::size_t y = 0; y < mat.rows(); ++y){
            std::fprintf(file, format.c_str(), mat.at(y, x));
        }
    }
}

template <typename T>
void writeVector(std::FILE* file, const Vector<T>& vec, const std::string& formatCode){
    mm_write_mtx_array_size(file, 1, static_cast<int>(vec.size()));

    const std::string format = formatCode + "\n";

    for(const auto& i : vec){
        std::fprintf(file, format.c_str(), i);
    }
}

int validFile(const std::string& fileName){
    std::FILE* file = std::fopen(fileName.c_str(), "r");

    //Unable to open
    if(!file) return 1;

    MM_typecode code{};
    if(mm_read_banner(file, &code)){
        //Error trying to read the banner
        std::fclose(file);
        return 2;
    }

    bool valid = true;

    //Ensure file contains a matrix
    if(!mm_is_matrix(code))     valid = false;

    //Ensure matrix is dense
    if(!mm_is_dense(code))      valid = false;

    //Ensure matrix is array format, not coordinate
    if(!mm_is_array(code))      valid = false;

    //Ensure matrix is general (not symmetric or something weird)
    if(!mm_is_general(code))    valid = false;

    //Ensure matrix either contains integers or real(floating point)
    if(!mm_is_real(code) && !mm_is_integer(code))   valid = false;

    std::fclose(file);

    if(valid) return 0;

    return 3;
}

DataType readDataType(const std::string& fileName){
    std::FILE* file = std::fopen(fileName.c_str(), "r");
    MM_typecode code{};
    mm_read_banner(file, &code);
    std::fclose(file);

    if(mm_is_real(code))    return DataType::real;
    if(mm_is_integer(code)) return DataType::integer;
    assert(false);
}

Matrix<int> readMatrixInt(std::FILE* file){
    const std::string format = "%i";
    return readMatrix<int>(file, format);
}

Matrix<double> readMatrixDouble(std::FILE* file){
    const std::string format = "%lg";
    return readMatrix<double>(file, format);
}

Vector<int> readVectorInt(std::FILE* file){
    const std::string format = "%i";
    return readVector<int>(file, format);
}

Vector<double> readVectorDouble(std::FILE* file){
    const std::string format = "%lg";
    return readVector<double>(file, format);
}


bool writeMatrixInt(const std::string& fileName, const Matrix<int>& mat){
    std::FILE* file = std::fopen(fileName.c_str(), "w");
    if(!file) return false;

    const std::string format = "%i";
    writeHeader(file, DataType::integer);
    writeMatrix(file, mat, format);

    std::fclose(file);
    return true;
}

bool writeMatrixDouble(const std::string& fileName, const Matrix<double>& mat){
    std::FILE* file = std::fopen(fileName.c_str(), "w");
    if(!file) return false;

    const std::string format = "%lg";
    writeHeader(file, DataType::real);
    writeMatrix(file, mat, format);

    std::fclose(file);
    return true;
}


bool writeVectorInt(const std::string& fileName, const Vector<int>& vec){
    std::FILE* file = std::fopen(fileName.c_str(), "w");
    if(!file) return false;

    const std::string format = "%i";
    writeHeader(file, DataType::integer);
    writeVector(file, vec, format);

    std::fclose(file);
    return true;
}

bool writeVectorDouble(const std::string& fileName, const Vector<double>& vec){
    std::FILE* file = std::fopen(fileName.c_str(), "w");
    if(!file) return false;

    const std::string format = "%lg";
    writeHeader(file, DataType::real);
    writeVector(file, vec, format);

    std::fclose(file);
    return true;
}

bool writeVTK(CellGrid<double, 2> grid, Matrix<double> matrix, std::string& file_name, std::string& grid_name, std::string& data_name) {

    std::FILE* file = std::fopen(file_name.c_str(), "w");
    if (!file) return false;

    // Write VTK header
    fprintf(file, "# vtk DataFile Version 2.0\n");

    // Write grid data
    fprintf(file, "%s\n", grid_name.c_str());
    fprintf(file, "ASCII\n");

    fprintf(file, "DATASET STRUCTURED_POINTS\n");
    fprintf(file, "DIMENSIONS %i %i %i\n",
        grid.N_pts[X],
        grid.N_pts[Y],
        1
    );
    fprintf(file, "ORIGIN %f %f %f\n",
        grid.lower_limit[X] + grid.spacing[X]/2,
        grid.lower_limit[Y] + grid.spacing[Y]/2,
        0.
    );
    fprintf(file, "SPACING %f %f %f\n",
        grid.spacing[X],
        grid.spacing[Y],
        1.
    );

    // Write data
    fprintf(file, "\nPOINT_DATA %i\n", grid.N_pts[X]*grid.N_pts[Y]);
    fprintf(file, "SCALARS %s double\n", data_name.c_str());
    fprintf(file, "LOOKUP_TABLE default\n");
    for (int i = 0; i < grid.N_pts[X]; i++) {
        for (int j = 0; j < grid.N_pts[Y]; j++) {
            fprintf(file, "%f ", matrix.at(i,j));
        }
        fprintf(file, "\n");
    }

    fclose(file);
    return true;

}


} // END NAMESPACE hps
