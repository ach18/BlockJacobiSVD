#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
#include "src/utils/types.hpp"
#include "src/utils/util.hpp"
#include "src/svd/two-sided/svd.hpp"

int main(int argc, char* argv[]) 
{
    size_t n;
    size_t block_size;

    if (argc < 3)
    {
        std::cout << "First argument must be a size of matrix" << std::endl;
        std::cout << "Second argument must be a block size" << std::endl;
        return 1;
    }

    n = atoi(argv[1]);
    block_size = atoi(argv[2]);
    std::cout << "Singular decomposition of double square matrix" << std::endl;

    aligned_vector<double> A(n * n);
    aligned_vector<double> B(n * n);
    aligned_vector<double> U(n * n, 0);
    aligned_vector<double> V(n * n, 0);

    matrix_t Data_matr = {&A[0], n, n};
    matrix_t B_mat = {&B[0], n, n};
    matrix_t U_mat = {&U[0], n, n};
    matrix_t V_mat = {&V[0], n, n};

    matrix_from_file(Data_matr, "./file.in");
    assert(Data_matr.rows == Data_matr.cols);
    assert(Data_matr.rows == B_mat.rows && Data_matr.cols == B_mat.cols);
    assert(Data_matr.rows == U_mat.rows && Data_matr.cols == U_mat.cols);
    assert(Data_matr.rows == V_mat.rows && Data_matr.cols == V_mat.cols);

    size_t iterations = svd_blocked(Data_matr, B_mat, U_mat, V_mat, block_size);

    matrix_to_file(B_mat, "./AFile.to");
    matrix_to_file(U_mat, "./UFile.to");
    matrix_to_file(V_mat, "./VFile.to");
    
    return 0;
}