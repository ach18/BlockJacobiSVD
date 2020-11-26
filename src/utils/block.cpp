#include "block.hpp"
#include <immintrin.h>

void copy_block(struct matrix_t Smat, size_t blockS_row, size_t blockS_col, struct matrix_t Dmat, size_t blockD_row,
                size_t blockD_col, size_t block_size) {
    size_t nS = Smat.rows;
    size_t nD = Dmat.rows;
    double* S = Smat.ptr;
    double* D = Dmat.ptr;
    size_t Sbeg = blockS_row * block_size * nS + blockS_col * block_size;
    size_t Dbeg = blockD_row * block_size * nD + blockD_col * block_size;
    for (size_t i = 0; i < block_size; ++i) {
        for (size_t j = 0; j < block_size; j += 8) {
            D[Dbeg + i * nD + j + 0] = S[Sbeg + i * nS + j + 0];
            D[Dbeg + i * nD + j + 1] = S[Sbeg + i * nS + j + 1];
            D[Dbeg + i * nD + j + 2] = S[Sbeg + i * nS + j + 2];
            D[Dbeg + i * nD + j + 3] = S[Sbeg + i * nS + j + 3];
            D[Dbeg + i * nD + j + 4] = S[Sbeg + i * nS + j + 4];
            D[Dbeg + i * nD + j + 5] = S[Sbeg + i * nS + j + 5];
            D[Dbeg + i * nD + j + 6] = S[Sbeg + i * nS + j + 6];
            D[Dbeg + i * nD + j + 7] = S[Sbeg + i * nS + j + 7];
        }
    }
}

static inline void inner_product(size_t n, double* A, size_t lda, double* B, size_t ldb, double* C, size_t ldc) {
    double* C_1 = C;
    double* C_2 = C_1 + ldc;
    double* C_3 = C_2 + ldc;
    double* C_4 = C_3 + ldc;

    double* A_1 = A;
    double* A_2 = A_1 + lda;
    double* A_3 = A_2 + lda;
    double* A_4 = A_3 + lda;

    __m256d c[8];

    c[0] = _mm256_setzero_pd();
    c[1] = _mm256_setzero_pd();
    c[2] = _mm256_setzero_pd();
    c[3] = _mm256_setzero_pd();
    c[4] = _mm256_setzero_pd();
    c[5] = _mm256_setzero_pd();
    c[6] = _mm256_setzero_pd();
    c[7] = _mm256_setzero_pd();

    for (size_t p = 0; p < n; p += 1) {
        __m256d b = _mm256_load_pd(B);
        __m256d b1 = _mm256_load_pd(B + 4);
        B = B + ldb;

        __m256d a1 = _mm256_broadcast_sd(A_1++);
        __m256d a2 = _mm256_broadcast_sd(A_2++);
        __m256d a3 = _mm256_broadcast_sd(A_3++);
        __m256d a4 = _mm256_broadcast_sd(A_4++);

        c[0] = _mm256_fmadd_pd(a1, b, c[0]);
        c[1] = _mm256_fmadd_pd(a1, b1, c[1]);
        c[2] = _mm256_fmadd_pd(a2, b, c[2]);
        c[3] = _mm256_fmadd_pd(a2, b1, c[3]);
        c[4] = _mm256_fmadd_pd(a3, b, c[4]);
        c[5] = _mm256_fmadd_pd(a3, b1, c[5]);
        c[6] = _mm256_fmadd_pd(a4, b, c[6]);
        c[7] = _mm256_fmadd_pd(a4, b1, c[7]);
    }

    _mm256_store_pd(C_1, c[0]);
    _mm256_store_pd(C_1 + 4, c[1]);
    _mm256_store_pd(C_2, c[2]);
    _mm256_store_pd(C_2 + 4, c[3]);
    _mm256_store_pd(C_3, c[4]);
    _mm256_store_pd(C_3 + 4, c[5]);
    _mm256_store_pd(C_4, c[6]);
    _mm256_store_pd(C_4 + 4, c[7]);
}

static inline void inner_add_product(size_t n, double* A, size_t lda, double* B, size_t ldb, double* C, size_t ldc,
                                     double* D, size_t ldd) {
    double* A_1 = A;
    double* A_2 = A_1 + lda;
    double* A_3 = A_2 + lda;
    double* A_4 = A_3 + lda;

    double* C_1 = C;
    double* C_2 = C_1 + ldc;
    double* C_3 = C_2 + ldc;
    double* C_4 = C_3 + ldc;

    double* D_1 = D;
    double* D_2 = D_1 + ldd;
    double* D_3 = D_2 + ldd;
    double* D_4 = D_3 + ldd;

    __m256d c[8];

    c[0] = _mm256_load_pd(C_1);
    c[1] = _mm256_load_pd(C_1 + 4);
    c[2] = _mm256_load_pd(C_2);
    c[3] = _mm256_load_pd(C_2 + 4);
    c[4] = _mm256_load_pd(C_3);
    c[5] = _mm256_load_pd(C_3 + 4);
    c[6] = _mm256_load_pd(C_4);
    c[7] = _mm256_load_pd(C_4 + 4);

    for (size_t p = 0; p < n; p += 1) {
        __m256d b = _mm256_load_pd(B);
        __m256d b1 = _mm256_load_pd(B + 4);
        B = B + ldb;

        __m256d a1 = _mm256_broadcast_sd(A_1++);
        __m256d a2 = _mm256_broadcast_sd(A_2++);
        __m256d a3 = _mm256_broadcast_sd(A_3++);
        __m256d a4 = _mm256_broadcast_sd(A_4++);

        c[0] = _mm256_fmadd_pd(a1, b, c[0]);
        c[1] = _mm256_fmadd_pd(a1, b1, c[1]);
        c[2] = _mm256_fmadd_pd(a2, b, c[2]);
        c[3] = _mm256_fmadd_pd(a2, b1, c[3]);
        c[4] = _mm256_fmadd_pd(a3, b, c[4]);
        c[5] = _mm256_fmadd_pd(a3, b1, c[5]);
        c[6] = _mm256_fmadd_pd(a4, b, c[6]);
        c[7] = _mm256_fmadd_pd(a4, b1, c[7]);
    }

    _mm256_store_pd(D_1, c[0]);
    _mm256_store_pd(D_1 + 4, c[1]);
    _mm256_store_pd(D_2, c[2]);
    _mm256_store_pd(D_2 + 4, c[3]);
    _mm256_store_pd(D_3, c[4]);
    _mm256_store_pd(D_3 + 4, c[5]);
    _mm256_store_pd(D_4, c[6]);
    _mm256_store_pd(D_4 + 4, c[7]);
}

// perform C = AB
void mult_block(struct matrix_t Amat, size_t blockA_row, size_t blockA_col, struct matrix_t Bmat, size_t blockB_row,
                size_t blockB_col, struct matrix_t Cmat, size_t blockC_row, size_t blockC_col, size_t block_size) {
    size_t nA = Amat.rows;
    size_t nB = Bmat.rows;
    size_t nC = Cmat.rows;
    double* A = Amat.ptr;
    double* B = Bmat.ptr;
    double* C = Cmat.ptr;
    size_t Abeg = blockA_row * block_size * nA + blockA_col * block_size;
    size_t Bbeg = blockB_row * block_size * nB + blockB_col * block_size;
    size_t Cbeg = blockC_row * block_size * nC + blockC_col * block_size;

    for (size_t i = 0; i < block_size; i += 8) {
        for (size_t j = 0; j < block_size; j += 4) {
            inner_product(block_size, &A[Abeg + j * nA], nA, &B[Bbeg + i], nB, &C[Cbeg + j * nC + i], nC);
        }
    }
}

// perform D = C + AB
void mult_add_block(struct matrix_t Amat, size_t blockA_row, size_t blockA_col, struct matrix_t Bmat, size_t blockB_row,
                    size_t blockB_col, struct matrix_t Cmat, size_t blockC_row, size_t blockC_col, struct matrix_t Dmat,
                    size_t blockD_row, size_t blockD_col, size_t block_size) {
    size_t nA = Amat.rows;
    size_t nB = Bmat.rows;
    size_t nC = Cmat.rows;
    size_t nD = Dmat.rows;
    double* A = Amat.ptr;
    double* B = Bmat.ptr;
    double* C = Cmat.ptr;
    double* D = Dmat.ptr;
    size_t Abeg = blockA_row * block_size * nA + blockA_col * block_size;
    size_t Bbeg = blockB_row * block_size * nB + blockB_col * block_size;
    size_t Cbeg = blockC_row * block_size * nC + blockC_col * block_size;
    size_t Dbeg = blockD_row * block_size * nD + blockD_col * block_size;

    for (size_t i = 0; i < block_size; i += 8) {
        for (size_t j = 0; j < block_size; j += 4) {
            inner_add_product(block_size, &A[Abeg + j * nA], nA, &B[Bbeg + i], nB, &C[Cbeg + j * nC + i], nC,
                              &D[Dbeg + j * nD + i], nD);
        }
    }
}

static inline void inner_product_left_tr(size_t n, double* A, size_t lda, double* B, size_t ldb, double* C,
                                         size_t ldc) {
    double* C_1 = C;
    double* C_2 = C_1 + ldc;
    double* C_3 = C_2 + ldc;
    double* C_4 = C_3 + ldc;

    __m256d c[8];

    c[0] = _mm256_setzero_pd();
    c[1] = _mm256_setzero_pd();
    c[2] = _mm256_setzero_pd();
    c[3] = _mm256_setzero_pd();
    c[4] = _mm256_setzero_pd();
    c[5] = _mm256_setzero_pd();
    c[6] = _mm256_setzero_pd();
    c[7] = _mm256_setzero_pd();

    for (size_t p = 0; p < n; p += 1) {
        __m256d b = _mm256_load_pd(B);
        __m256d b1 = _mm256_load_pd(B + 4);
        B = B + ldb;

        __m256d a = _mm256_load_pd(A);
        __m256d a1 = _mm256_permute4x64_pd(a, 0);
        __m256d a2 = _mm256_permute4x64_pd(a, 85);
        __m256d a3 = _mm256_permute4x64_pd(a, 170);
        __m256d a4 = _mm256_permute4x64_pd(a, 255);
        A = A + lda;

        c[0] = _mm256_fmadd_pd(a1, b, c[0]);
        c[1] = _mm256_fmadd_pd(a1, b1, c[1]);
        c[2] = _mm256_fmadd_pd(a2, b, c[2]);
        c[3] = _mm256_fmadd_pd(a2, b1, c[3]);
        c[4] = _mm256_fmadd_pd(a3, b, c[4]);
        c[5] = _mm256_fmadd_pd(a3, b1, c[5]);
        c[6] = _mm256_fmadd_pd(a4, b, c[6]);
        c[7] = _mm256_fmadd_pd(a4, b1, c[7]);
    }

    _mm256_store_pd(C_1, c[0]);
    _mm256_store_pd(C_1 + 4, c[1]);
    _mm256_store_pd(C_2, c[2]);
    _mm256_store_pd(C_2 + 4, c[3]);
    _mm256_store_pd(C_3, c[4]);
    _mm256_store_pd(C_3 + 4, c[5]);
    _mm256_store_pd(C_4, c[6]);
    _mm256_store_pd(C_4 + 4, c[7]);
}

static inline void inner_add_product_left_tr(size_t n, double* A, size_t lda, double* B, size_t ldb, double* C,
                                             size_t ldc, double* D, size_t ldd) {
    double* C_1 = C;
    double* C_2 = C_1 + ldc;
    double* C_3 = C_2 + ldc;
    double* C_4 = C_3 + ldc;

    double* D_1 = D;
    double* D_2 = D_1 + ldd;
    double* D_3 = D_2 + ldd;
    double* D_4 = D_3 + ldd;

    __m256d c[8];

    c[0] = _mm256_load_pd(C_1);
    c[1] = _mm256_load_pd(C_1 + 4);
    c[2] = _mm256_load_pd(C_2);
    c[3] = _mm256_load_pd(C_2 + 4);
    c[4] = _mm256_load_pd(C_3);
    c[5] = _mm256_load_pd(C_3 + 4);
    c[6] = _mm256_load_pd(C_4);
    c[7] = _mm256_load_pd(C_4 + 4);

    for (size_t p = 0; p < n; p += 1) {
        __m256d b = _mm256_load_pd(B);
        __m256d b1 = _mm256_load_pd(B + 4);
        B = B + ldb;

        __m256d a = _mm256_load_pd(A);
        __m256d a1 = _mm256_permute4x64_pd(a, 0);
        __m256d a2 = _mm256_permute4x64_pd(a, 85);
        __m256d a3 = _mm256_permute4x64_pd(a, 170);
        __m256d a4 = _mm256_permute4x64_pd(a, 255);
        A = A + lda;

        c[0] = _mm256_fmadd_pd(a1, b, c[0]);
        c[1] = _mm256_fmadd_pd(a1, b1, c[1]);
        c[2] = _mm256_fmadd_pd(a2, b, c[2]);
        c[3] = _mm256_fmadd_pd(a2, b1, c[3]);
        c[4] = _mm256_fmadd_pd(a3, b, c[4]);
        c[5] = _mm256_fmadd_pd(a3, b1, c[5]);
        c[6] = _mm256_fmadd_pd(a4, b, c[6]);
        c[7] = _mm256_fmadd_pd(a4, b1, c[7]);
    }

    _mm256_store_pd(D_1, c[0]);
    _mm256_store_pd(D_1 + 4, c[1]);
    _mm256_store_pd(D_2, c[2]);
    _mm256_store_pd(D_2 + 4, c[3]);
    _mm256_store_pd(D_3, c[4]);
    _mm256_store_pd(D_3 + 4, c[5]);
    _mm256_store_pd(D_4, c[6]);
    _mm256_store_pd(D_4 + 4, c[7]);
}

// perform C = (A^T)B
// for C_ij, use ith column of A and jth column of B
void mult_transpose_block(struct matrix_t Amat, size_t blockA_row, size_t blockA_col, struct matrix_t Bmat,
                          size_t blockB_row, size_t blockB_col, struct matrix_t Cmat, size_t blockC_row,
                          size_t blockC_col, size_t block_size) {
    size_t nA = Amat.rows;
    size_t nB = Bmat.rows;
    size_t nC = Cmat.rows;
    double* A = Amat.ptr;
    double* B = Bmat.ptr;
    double* C = Cmat.ptr;
    size_t Abeg = blockA_row * block_size * nA + blockA_col * block_size;
    size_t Bbeg = blockB_row * block_size * nB + blockB_col * block_size;
    size_t Cbeg = blockC_row * block_size * nC + blockC_col * block_size;

    for (size_t i = 0; i < block_size; i += 8) {
        for (size_t j = 0; j < block_size; j += 4) {
            inner_product_left_tr(block_size, &A[Abeg + j], nA, &B[Bbeg + i], nB, &C[Cbeg + j * nC + i], nC);
        }
    }
}

// perform D = C + (A^T)B
// for D_ij, use ith column of A and jth column of B.
void mult_add_transpose_block(struct matrix_t Amat, size_t blockA_row, size_t blockA_col, struct matrix_t Bmat,
                              size_t blockB_row, size_t blockB_col, struct matrix_t Cmat, size_t blockC_row,
                              size_t blockC_col, struct matrix_t Dmat, size_t blockD_row, size_t blockD_col,
                              size_t block_size) {
    size_t nA = Amat.rows;
    size_t nB = Bmat.rows;
    size_t nC = Cmat.rows;
    size_t nD = Dmat.rows;
    double* A = Amat.ptr;
    double* B = Bmat.ptr;
    double* C = Cmat.ptr;
    double* D = Dmat.ptr;
    size_t Abeg = blockA_row * block_size * nA + blockA_col * block_size;
    size_t Bbeg = blockB_row * block_size * nB + blockB_col * block_size;
    size_t Cbeg = blockC_row * block_size * nC + blockC_col * block_size;
    size_t Dbeg = blockD_row * block_size * nD + blockD_col * block_size;

    for (size_t i = 0; i < block_size; i += 8) {
        for (size_t j = 0; j < block_size; j += 4) {
            inner_add_product_left_tr(block_size, &A[Abeg + j], nA, &B[Bbeg + i], nB, &C[Cbeg + j * nC + i], nC,
                                      &D[Dbeg + j * nD + i], nD);
        }
    }
}

static inline __m128d hsum_double_avx(__m256d v) __attribute__((always_inline));
static inline __m128d hsum_double_avx(__m256d v) {
    __m128d vlow = _mm256_castpd256_pd128(v);
    __m128d vhigh = _mm256_extractf128_pd(v, 1);  // high 128
    vlow = _mm_add_pd(vlow, vhigh);               // reduce down to 128

    __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
    return _mm_add_pd(vlow, high64);  // reduce to scalar
}

static inline void inner_product_tr(size_t n, double* A, size_t lda, double* B, size_t ldb, double* C, size_t ldc) {
    double* A_1 = A;
    double* A_2 = A_1 + lda;
    double* A_3 = A_2 + lda;
    double* A_4 = A_3 + lda;

    double* B_1 = B;
    double* B_2 = B_1 + ldb;

    double* C_1 = C;
    double* C_2 = C_1 + ldc;
    double* C_3 = C_2 + ldc;
    double* C_4 = C_3 + ldc;

    __m256d c[8];

    c[0] = _mm256_setzero_pd();
    c[1] = _mm256_setzero_pd();
    c[2] = _mm256_setzero_pd();
    c[3] = _mm256_setzero_pd();
    c[4] = _mm256_setzero_pd();
    c[5] = _mm256_setzero_pd();
    c[6] = _mm256_setzero_pd();
    c[7] = _mm256_setzero_pd();

    for (size_t p = 0; p < n; p += 4) {
        __m256d a1 = _mm256_load_pd(A_1 + p);
        __m256d a2 = _mm256_load_pd(A_2 + p);
        __m256d a3 = _mm256_load_pd(A_3 + p);
        __m256d a4 = _mm256_load_pd(A_4 + p);

        __m256d b1 = _mm256_load_pd(B_1 + p);
        __m256d b2 = _mm256_load_pd(B_2 + p);

        c[0] = _mm256_fmadd_pd(a1, b1, c[0]);
        c[1] = _mm256_fmadd_pd(a1, b2, c[1]);
        c[2] = _mm256_fmadd_pd(a2, b1, c[2]);
        c[3] = _mm256_fmadd_pd(a2, b2, c[3]);
        c[4] = _mm256_fmadd_pd(a3, b1, c[4]);
        c[5] = _mm256_fmadd_pd(a3, b2, c[5]);
        c[6] = _mm256_fmadd_pd(a4, b1, c[6]);
        c[7] = _mm256_fmadd_pd(a4, b2, c[7]);
    }

    __m128d c00 = hsum_double_avx(c[0]);
    __m128d c01 = hsum_double_avx(c[1]);
    __m128d c10 = hsum_double_avx(c[2]);
    __m128d c11 = hsum_double_avx(c[3]);
    __m128d c20 = hsum_double_avx(c[4]);
    __m128d c21 = hsum_double_avx(c[5]);
    __m128d c30 = hsum_double_avx(c[6]);
    __m128d c31 = hsum_double_avx(c[7]);

    _mm_store_pd(C_1, _mm_shuffle_pd(c00, c01, 0));
    _mm_store_pd(C_2, _mm_shuffle_pd(c10, c11, 0));
    _mm_store_pd(C_3, _mm_shuffle_pd(c20, c21, 0));
    _mm_store_pd(C_4, _mm_shuffle_pd(c30, c31, 0));
}

static inline void inner_add_product_tr(size_t n, double* A, size_t lda, double* B, size_t ldb, double* C, size_t ldc,
                                        double* D, size_t ldd) {
    double* A_1 = A;
    double* A_2 = A_1 + lda;
    double* A_3 = A_2 + lda;
    double* A_4 = A_3 + lda;

    double* B_1 = B;
    double* B_2 = B_1 + ldb;

    double* C_1 = C;
    double* C_2 = C_1 + ldc;
    double* C_3 = C_2 + ldc;
    double* C_4 = C_3 + ldc;

    double* D_1 = D;
    double* D_2 = D_1 + ldd;
    double* D_3 = D_2 + ldd;
    double* D_4 = D_3 + ldd;

    __m256d c[8];

    c[0] = _mm256_setzero_pd();
    c[1] = _mm256_setzero_pd();
    c[2] = _mm256_setzero_pd();
    c[3] = _mm256_setzero_pd();
    c[4] = _mm256_setzero_pd();
    c[5] = _mm256_setzero_pd();
    c[6] = _mm256_setzero_pd();
    c[7] = _mm256_setzero_pd();

    for (size_t p = 0; p < n; p += 4) {
        __m256d a1 = _mm256_load_pd(A_1 + p);
        __m256d a2 = _mm256_load_pd(A_2 + p);
        __m256d a3 = _mm256_load_pd(A_3 + p);
        __m256d a4 = _mm256_load_pd(A_4 + p);

        __m256d b1 = _mm256_load_pd(B_1 + p);
        __m256d b2 = _mm256_load_pd(B_2 + p);

        c[0] = _mm256_fmadd_pd(a1, b1, c[0]);
        c[1] = _mm256_fmadd_pd(a1, b2, c[1]);
        c[2] = _mm256_fmadd_pd(a2, b1, c[2]);
        c[3] = _mm256_fmadd_pd(a2, b2, c[3]);
        c[4] = _mm256_fmadd_pd(a3, b1, c[4]);
        c[5] = _mm256_fmadd_pd(a3, b2, c[5]);
        c[6] = _mm256_fmadd_pd(a4, b1, c[6]);
        c[7] = _mm256_fmadd_pd(a4, b2, c[7]);
    }

    __m128d c00 = hsum_double_avx(c[0]);
    __m128d c01 = hsum_double_avx(c[1]);
    __m128d c10 = hsum_double_avx(c[2]);
    __m128d c11 = hsum_double_avx(c[3]);
    __m128d c20 = hsum_double_avx(c[4]);
    __m128d c21 = hsum_double_avx(c[5]);
    __m128d c30 = hsum_double_avx(c[6]);
    __m128d c31 = hsum_double_avx(c[7]);

    __m128d c0 = _mm_load_pd(C_1);
    __m128d c1 = _mm_load_pd(C_2);
    __m128d c2 = _mm_load_pd(C_3);
    __m128d c3 = _mm_load_pd(C_4);

    _mm_store_pd(D_1, _mm_add_pd(c0, _mm_shuffle_pd(c00, c01, 0)));
    _mm_store_pd(D_2, _mm_add_pd(c1, _mm_shuffle_pd(c10, c11, 0)));
    _mm_store_pd(D_3, _mm_add_pd(c2, _mm_shuffle_pd(c20, c21, 0)));
    _mm_store_pd(D_4, _mm_add_pd(c3, _mm_shuffle_pd(c30, c31, 0)));
}

// perform C = A(B^T)
void mult_transpose_block_right(struct matrix_t Amat, size_t blockA_row, size_t blockA_col, struct matrix_t Bmat,
                                size_t blockB_row, size_t blockB_col, struct matrix_t Cmat, size_t blockC_row,
                                size_t blockC_col, size_t block_size) {
    size_t nA = Amat.rows;
    size_t nB = Bmat.rows;
    size_t nC = Cmat.rows;
    double* A = Amat.ptr;
    double* B = Bmat.ptr;
    double* C = Cmat.ptr;
    size_t Abeg = blockA_row * block_size * nA + blockA_col * block_size;
    size_t Bbeg = blockB_row * block_size * nB + blockB_col * block_size;
    size_t Cbeg = blockC_row * block_size * nC + blockC_col * block_size;
    for (size_t i = 0; i < block_size; i += 4) {
        for (size_t j = 0; j < block_size; j += 2) {
            inner_product_tr(block_size, &A[Abeg + i * nA], nA, &B[Bbeg + j * nB], nB, &C[Cbeg + i * nC + j], nC);
        }
    }
}

// perform C = C + A(B^T)
void mult_add_transpose_block_right(struct matrix_t Amat, size_t blockA_row, size_t blockA_col, struct matrix_t Bmat,
                                    size_t blockB_row, size_t blockB_col, struct matrix_t Cmat, size_t blockC_row,
                                    size_t blockC_col, struct matrix_t Dmat, size_t blockD_row, size_t blockD_col,
                                    size_t block_size) {
    size_t nA = Amat.rows;
    size_t nB = Bmat.rows;
    size_t nC = Cmat.rows;
    size_t nD = Dmat.rows;
    double* A = Amat.ptr;
    double* B = Bmat.ptr;
    double* C = Cmat.ptr;
    double* D = Dmat.ptr;
    size_t Abeg = blockA_row * block_size * nA + blockA_col * block_size;
    size_t Bbeg = blockB_row * block_size * nB + blockB_col * block_size;
    size_t Cbeg = blockC_row * block_size * nC + blockC_col * block_size;
    size_t Dbeg = blockD_row * block_size * nD + blockD_col * block_size;
    for (size_t i = 0; i < block_size; i += 4) {
        for (size_t j = 0; j < block_size; j += 2) {
            inner_add_product_tr(block_size, &A[Abeg + i * nA], nA, &B[Bbeg + j * nB], nB, &C[Cbeg + i * nC + j], nC,
                                 &D[Dbeg + i * nD + j], nD);
        }
    }
}
