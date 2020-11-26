#pragma once
#include <cstdlib>
#include "matrix.hpp"

// perform C = AB
void mult_block(struct matrix_t Amat, size_t blockA_row, size_t blockA_col, struct matrix_t Bmat, size_t blockB_row,
                size_t blockB_col, struct matrix_t Cmat, size_t blockC_row, size_t blockC_col, size_t block_size);
// perform D = C + AB
void mult_add_block(struct matrix_t Amat, size_t blockA_row, size_t blockA_col, struct matrix_t Bmat, size_t blockB_row,
                    size_t blockB_col, struct matrix_t Cmat, size_t blockC_row, size_t blockC_col, struct matrix_t Dmat,
                    size_t blockD_row, size_t blockD_col, size_t block_size);
// perform C = (A^T)B
void mult_transpose_block(struct matrix_t Amat, size_t blockA_row, size_t blockA_col, struct matrix_t Bmat,
                          size_t blockB_row, size_t blockB_col, struct matrix_t Cmat, size_t blockC_row,
                          size_t blockC_col, size_t block_size);
// perform D = C + (A^T)B
void mult_add_transpose_block(struct matrix_t Amat, size_t blockA_row, size_t blockA_col, struct matrix_t Bmat,
                              size_t blockB_row, size_t blockB_col, struct matrix_t Cmat, size_t blockC_row,
                              size_t blockC_col, struct matrix_t Dmat, size_t blockD_row, size_t blockD_col,
                              size_t block_size);
// perform C = A(B^T)
void mult_transpose_block_right(struct matrix_t Amat, size_t blockA_row, size_t blockA_col, struct matrix_t Bmat,
                                size_t blockB_row, size_t blockB_col, struct matrix_t Cmat, size_t blockC_row,
                                size_t blockC_col, size_t block_size);
// perform D = C + A(B^T)
void mult_add_transpose_block_right(struct matrix_t Amat, size_t blockA_row, size_t blockA_col, struct matrix_t Bmat,
                                    size_t blockB_row, size_t blockB_col, struct matrix_t Cmat, size_t blockC_row,
                                    size_t blockC_col, struct matrix_t Dmat, size_t blockD_row, size_t blockD_col,
                                    size_t block_size);
void copy_block(struct matrix_t S, size_t blockS_row, size_t blockS_col, struct matrix_t D, size_t blockD_row,
                size_t blockD_col, size_t block_size);
