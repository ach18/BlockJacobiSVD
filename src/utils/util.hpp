#pragma once
#include <float.h>
#include "types.hpp"

/**
 * Compare two double precision floating point values and return true if they
 * are close in the relative scale.
 *
 * Formally, the function does true if |x - y| <= eps * |x + y| is satisfied.
 *
 * @param x Double value.
 * @param y Double value.
 * @param eps Epsilon value to check for closeness.
 */
bool isclose(double x, double y, double eps = DBL_EPSILON);

bool is_normalized(double v);

/**
 * Given three entries of a symmetric matrix X, compute a cosine-sine pair (c,
 * s)
 * such that if B = J(i,j,t)^T*A*J(i,j,t), then B_{ij} = B_{ji} = 0. Here,
 * J(i, j, t) represent the Jacobi rotation matrix constructed using (c, s).
 *
 * The algorithm is described in
 * http://web.mit.edu/ehliu/Public/sclark/Golub%20G.H.,%20Van%20Loan%20C.F.-%20Matrix%20Computations.pdf
 * as Algorithm 8.4.1.
 *
 * @param x_ii Top left entry of the block of X, i.e. X[i][i].
 * @param x_ij Cross entry of the block of X, i.e. X[i][j].
 * @param x_jj Bottom right entry of the block of X, i.e. X[j][j].
 * @param[out] c Cosine value.
 * @param[out] s Sine value.
 */
void sym_jacobi_coeffs(double x_ii, double x_ij, double x_jj, double* c, double* s);

/**
 * Comparator function that returns 1 (true) if the first element comes before
 * the second element in the ordering implied by the comparison.
 *
 * For example, if the function implements '<' operator, and the inputs are
 * 3 and 5, the function must return 1. In the same case, if the inputs are
 * 8 and 3, the function must return 0. Similary, if the inputs are 5 and 5,
 * the function must return 0.
 */
typedef int (*comparator)(double, double);

/**
 * Comparator that returns 1 if (x < y), and 0 otherwise.
 */
int less(double x, double y);

/**
 * Comparator that returns 1 if (x > y), and 0 otherwise.
 */
int greater(double x, double y);

/**
 * Return 1 if the number if positive and -1 otherwise.
 */
double sign(double x);

/**
 * Reorder a given matrix decomposition according to the singular/eigen values
 * and a comparison function.
 *
 * @param vals Singular/eigen values of the decomposition.
 * @param An array of matrix types. For each matrix in the array, ith column of
 * that matrix must correspond to the ith singular/eigen value. These columns
 * will be reordered in the same manner as their corresponding singular/eigen
 * values.
 * @param n_matrices Number of matrices in the given matrix array.
 * @param cmp_fn A comparator to use when reordering the columns. This function
 * must be a total ordering, and must define the sorting order. For example,
 * if cmp_fn implements less-than operator, then the values will be ordered in
 * ascending order.
 */
void reorder_decomposition(struct vector_t vals, struct matrix_t* matrices, int n_matrices, comparator cmp_fn);

void matrix_from_file(matrix_t A, char path[]);
void matrix_to_file(matrix_t A, char path[]);