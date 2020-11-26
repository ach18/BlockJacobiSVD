#pragma once

#include "types.hpp"

/**
 * @brief Initialize an identity matrix
 * @param Pmat Input matrix
 */
void matrix_identity(matrix_t Pmat);

/**
 * @brief Generate the transpose of a matrix
 * @param Pmat Input matrix
 * @param Qmat Output matrix - Transpose of input matrix
 */
void matrix_transpose(matrix_t Pmat, matrix_t Qmat);

/**
 * @brief Multiply two matrices and store the result in another matrix
 * @param Qmat Input left matrix
 * @param Rmat Input right matrix
 * @param Pmat Output matrix
 */
void matrix_mult(matrix_t Pmat, matrix_t Qmat, matrix_t Rmat);

// C = A + B
void matrix_add(matrix_t Amat, matrix_t Bmat, matrix_t Cmat);

/**
 * @brief Copies a matrix Amat to matrix Bmat
 * @param Bmat Output matrix, copy of input matrix
 * @param Amat Input matrix
 */

void matrix_copy(matrix_t Bmat, matrix_t Amat);

/**
 * @brief Calcluates the oof frebonius norm of matrix m
 * @param m Intput matrix
 * @param off_norm, output value, the off diagonal frbonius norm.
 */

void matrix_off_frobenius(matrix_t m, double* off_norm);

/**
 * @brief Calcluates the frebonius norm of matrix m
 * @param m Intput matrix
 * @param off_norm, output value, the frbonius norm.
 */

void matrix_frobenius(matrix_t m, double* norm, double* off_norm);
