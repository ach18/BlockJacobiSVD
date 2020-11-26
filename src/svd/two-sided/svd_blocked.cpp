#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include "../../utils/block.hpp"
#include "../../utils/matrix.hpp"
#include "nsvd.hpp"
#include "svd.hpp"
#include "svd_subprocedure.hpp"
#include "../../utils/types.hpp"
#include "../../utils/util.hpp"

/**
 * @param matrix_t Amat симметричная матрица A
 * @param matrix_t Bmat матрица A после деления на блоки
 * @param matrix_t Umat матрица левых сингулярных векторов
 * @param matrix_t Vmat матрица правых сингулярных векторов
 * @param size_t block_size размер блока матрицы Bmat
 * @return size_t sweeps число разверток методом Якоби
 **/
size_t svd_blocked(struct matrix_t Amat, struct matrix_t Bmat, struct matrix_t Umat, struct matrix_t Vmat,
                   size_t block_size) {

    size_t sweeps = 0;  //число повторений цикла развертки
    size_t block_iter = 0;
    const double tol = 1e-10;  //точность предела сходимости
    const size_t n = Amat.rows; //размер матрицы
    double norm = 0.0;      //норма Фробениуса матрицы B
    double off_norm = 0.0;  //норма Фробениуса только недиагональных элементов матрицы B

    matrix_copy(Bmat, Amat); //инициализация рабочей матрицы B
    //матрицы U и V инициализируются как единичные матрицы
    matrix_identity(Umat);
    matrix_identity(Vmat);
    matrix_frobenius(Bmat, &norm, &off_norm);

    const size_t n_blocks = n / block_size; //число блоков вдоль строки/столбца

    //если строка/столбец состоят из двух блоков, лучше вычислить SVD классическим не блочным Якоби
    if (n <= 2 * block_size) {
        size_t block_iters = svd_subprocedure(Bmat, Umat, Vmat);
        return block_iters;
    }

    //общее число элементов всех блоков строки/столбца должно быть равно размерности матрицы
    assert(n_blocks * block_size == n);

    //выделение памяти для хранения блоков матриц B U V M1 M2
    double* memory_block = (double*) malloc((4 + 4 + 4 + 1 + 1) * block_size * block_size * sizeof(double));
    double* Bblock = memory_block;
    double* Ublock = Bblock + 4 * block_size * block_size;
    double* Vblock = Ublock + 4 * block_size * block_size;
    //M1 M2 хранят промежуточные значения вычислений
    double* M1 = Vblock + 4 * block_size * block_size;
    double* M2 = M1 + block_size * block_size;

    //хранение блоков в виде структур matrix_t
    matrix_t Bblockmat = {Bblock, 2 * block_size, 2 * block_size};
    matrix_t Ublockmat = {Ublock, 2 * block_size, 2 * block_size};
    matrix_t Vblockmat = {Vblock, 2 * block_size, 2 * block_size};
    matrix_t M1mat = {M1, block_size, block_size};
    matrix_t M2mat = {M2, block_size, block_size};

    //основной цикл развретки, он продолжается пока 
    while (sqrt(off_norm) > tol * sqrt(norm)) {
        //цикл обхода над/поддиагональных блоков
        for (size_t i_block = 0; i_block < n_blocks - 1; ++i_block) {
            for (size_t j_block = i_block + 1; j_block < n_blocks; ++j_block) {
                //в Bblockmat копируются блоки с индексами ii ij ji jj из матрицы Bmat
                copy_block(Bmat, i_block, i_block, Bblockmat, 0, 0, block_size);
                copy_block(Bmat, i_block, j_block, Bblockmat, 0, 1, block_size);
                copy_block(Bmat, j_block, i_block, Bblockmat, 1, 0, block_size);
                copy_block(Bmat, j_block, j_block, Bblockmat, 1, 1, block_size);

                //вычисление SVD разложения циклическим методом Якоби над блоком Bblockmat
                block_iter += svd_subprocedure(Bblockmat, Ublockmat, Vblockmat);

                //транспонировать блок матрицы U 
                matrix_transpose(Ublockmat, Ublockmat);

                //обновление блоков матрицы B по строкам i j
                for (size_t k_block = 0; k_block < n_blocks; ++k_block) {
                    mult_block(Ublockmat, 0, 0, Bmat, i_block, k_block, M1mat, 0, 0, block_size);
                    mult_block(Ublockmat, 0, 1, Bmat, j_block, k_block, M2mat, 0, 0, block_size);
                    matrix_add(M1mat, M2mat, M2mat);
                    mult_block(Ublockmat, 1, 0, Bmat, i_block, k_block, M1mat, 0, 0, block_size);
                    copy_block(M2mat, 0, 0, Bmat, i_block, k_block, block_size);
                    mult_block(Ublockmat, 1, 1, Bmat, j_block, k_block, M2mat, 0, 0, block_size);
                    matrix_add(M1mat, M2mat, M2mat);
                    copy_block(M2mat, 0, 0, Bmat, j_block, k_block, block_size);
                }

                //обновление блоков матрицы B по столбцам i j
                for (size_t k_block = 0; k_block < n_blocks; ++k_block) {
                    mult_block(Bmat, k_block, i_block, Vblockmat, 0, 0, M1mat, 0, 0, block_size);
                    mult_block(Bmat, k_block, j_block, Vblockmat, 1, 0, M2mat, 0, 0, block_size);
                    matrix_add(M1mat, M2mat, M2mat);
                    mult_block(Bmat, k_block, i_block, Vblockmat, 0, 1, M1mat, 0, 0, block_size);
                    copy_block(M2mat, 0, 0, Bmat, k_block, i_block, block_size);
                    mult_block(Bmat, k_block, j_block, Vblockmat, 1, 1, M2mat, 0, 0, block_size);
                    matrix_add(M1mat, M2mat, M2mat);
                    copy_block(M2mat, 0, 0, Bmat, k_block, j_block, block_size);
                }

                //возвращение блока матрицы U^T в исходное состояние U 
                matrix_transpose(Ublockmat, Ublockmat);

                //обновление блоков матрицы U по столбцам i j
                for (size_t k_block = 0; k_block < n_blocks; ++k_block) {
                    mult_block(Umat, k_block, i_block, Ublockmat, 0, 0, M1mat, 0, 0, block_size);
                    mult_block(Umat, k_block, j_block, Ublockmat, 1, 0, M2mat, 0, 0, block_size);
                    matrix_add(M1mat, M2mat, M2mat);
                    mult_block(Umat, k_block, i_block, Ublockmat, 0, 1, M1mat, 0, 0, block_size);
                    copy_block(M2mat, 0, 0, Umat, k_block, i_block, block_size);
                    mult_block(Umat, k_block, j_block, Ublockmat, 1, 1, M2mat, 0, 0, block_size);
                    matrix_add(M1mat, M2mat, M2mat);
                    copy_block(M2mat, 0, 0, Umat, k_block, j_block, block_size);
                }

                //обновление блоков матрицы V по столбцам i j
                for (size_t k_block = 0; k_block < n_blocks; ++k_block) {
                    mult_block(Vmat, k_block, i_block, Vblockmat, 0, 0, M1mat, 0, 0, block_size);
                    mult_block(Vmat, k_block, j_block, Vblockmat, 1, 0, M2mat, 0, 0, block_size);
                    matrix_add(M1mat, M2mat, M2mat);
                    mult_block(Vmat, k_block, i_block, Vblockmat, 0, 1, M1mat, 0, 0, block_size);
                    copy_block(M2mat, 0, 0, Vmat, k_block, i_block, block_size);
                    mult_block(Vmat, k_block, j_block, Vblockmat, 1, 1, M2mat, 0, 0, block_size);
                    matrix_add(M1mat, M2mat, M2mat);
                    copy_block(M2mat, 0, 0, Vmat, k_block, j_block, block_size);
                }
            }
        }

        matrix_frobenius(Bmat, &norm, &off_norm);
        sweeps++;
    }

    free(memory_block);

    return sweeps;
}
