#include "matrix.hpp"
#include <immintrin.h>
#include <algorithm>
#include <cassert>
#include <utility>

void matrix_identity(matrix_t Pmat) {
    double* P = Pmat.ptr;
    const size_t M = Pmat.rows;
    const size_t N = Pmat.cols;

    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < N; j++) {
            if (i == j) {
                P[i * N + j] = 1.0;
            } else {
                P[i * N + j] = 0.0;
            }
        }
    }
}

void matrix_transpose(matrix_t Pmat, matrix_t Qmat) {
    double* P = Pmat.ptr;
    double* Q = Qmat.ptr;
    const size_t M = Pmat.rows;
    const size_t N = Pmat.cols;

    // Swap for in-space transpose
    if (P == Q) {
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < i; j++) {
                std::swap(Q[i * N + j], Q[j * N + i]);
            }
        }
    } else {
        for (size_t i = 0; i < M; i++) {
            for (size_t j = 0; j < N; j++) {
                Q[j * N + i] = P[i * N + j];
            }
        }
    }
}

void matrix_mult(matrix_t Pmat, matrix_t Qmat, matrix_t Rmat) {
    double* P = Pmat.ptr;
    double* Q = Qmat.ptr;
    double* R = Rmat.ptr;
    const size_t M = Pmat.rows;
    const size_t N = Pmat.cols;
    const size_t O = Qmat.cols;

    assert(Qmat.cols == Rmat.rows);

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            double sum = 0.0;

            for (size_t k = 0; k < O; k++) {
                sum += Q[i * N + k] * R[k * N + j];
            }

            P[i * N + j] = sum;
        }
    }
}

void matrix_add(matrix_t Amat, matrix_t Bmat, matrix_t Cmat) {
    size_t n = Amat.rows;
    double* A = Amat.ptr;
    double* B = Bmat.ptr;
    double* C = Cmat.ptr;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            C[n * i + j] = A[n * i + j] + B[n * i + j];
        }
    }
}

void matrix_copy(matrix_t Bmat, matrix_t Amat) {
    assert(Amat.rows == Bmat.rows);
    assert(Amat.cols == Amat.cols);

    std::copy(Amat.ptr, Amat.ptr + Amat.rows * Amat.cols, Bmat.ptr);
}

void matrix_frobenius(matrix_t m, double* norm, double* off_norm) {
    const size_t M = m.rows;
    const size_t N = m.cols;
    double* data = m.ptr;
    double elems_sum = 0.0;           // sum m[i][j]^2 for 0 < i < M and 0 < j < N
    double off_diag_elems_sum = 0.0;  // sum m[i][j]^2 for 0 < i < M and 0 < j < N and i == j

    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            elems_sum += data[N * i + j] * data[N * i + j];

            if (i == j) {
                continue;
            } else {
                off_diag_elems_sum += data[N * i + j] * data[N * i + j];
            }
        }
    }

    *norm = elems_sum;
    *off_norm = off_diag_elems_sum;
}

void matrix_off_frobenius(matrix_t m, double* off_norm) {
    const size_t M = m.rows;
    const size_t N = m.cols;
    double* data = m.ptr;
    double off_diag_elems_sum = 0.0;  // sum m[i][j]^2 for 0 < i < M and 0 < j < N and i == j

    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (i == j) {
                continue;
            } else {
                off_diag_elems_sum += data[N * i + j] * data[N * i + j];
            }
        }
    }

    *off_norm = off_diag_elems_sum;
}