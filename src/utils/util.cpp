#include "util.hpp"
#include <assert.h>
#include <math.h>
#include <limits>
#include "types.hpp"
#include <iostream>
#include <fstream>

bool isclose(double x, double y, double eps) { return fabs(x - y) <= eps * fabs(x + y); }

bool is_normalized(double v) { return fabs(v) >= std::numeric_limits<double>::min(); }

void sym_jacobi_coeffs(double x_ii, double x_ij, double x_jj, double* c, double* s) {
    if (!isclose(x_ij, 0)) {
        double tau, t, out_c;
        tau = (x_jj - x_ii) / (2 * x_ij);
        if (tau >= 0) {
            t = 1.0 / (tau + sqrt(1 + tau * tau));
        } else {
            t = -1.0 / (-tau + sqrt(1 + tau * tau));
        }
        out_c = 1.0 / sqrt(1 + t * t);
        *c = out_c;
        *s = t * out_c;
    } else {
        *c = 1.0;
        *s = 0.0;
    }
}

int less(double x, double y) { return x < y; }

int greater(double x, double y) { return x > y; }

double sign(double x) { return (x > 0.0) ? 1.0 : -1.0; }

void reorder_decomposition(struct vector_t vals, struct matrix_t* matrices, int n_matrices, comparator cmp_fn) {
    double* s = vals.ptr;
    const int n_vals = vals.len;
    for (int i = 0; i < n_vals; ++i) {
        double s_last = s[i];
        int i_last = i;
        for (int j = i + 1; j < n_vals; ++j) {
            if (cmp_fn(s[j], s_last)) {
                s_last = s[j];
                i_last = j;
            }
        }
        if (i_last != i) {
            double tmp;
            tmp = s[i];
            s[i] = s[i_last];
            s[i_last] = tmp;

            for (int j = 0; j < n_matrices; ++j) {
                int rows = matrices[j].rows;
                int cols = matrices[j].cols;
                double* M = matrices[j].ptr;
                for (int k = 0; k < rows; ++k) {
                    tmp = M[k * cols + i];
                    M[k * cols + i] = M[k * cols + i_last];
                    M[k * cols + i_last] = tmp;
                }
            }
        }
    }
}

void matrix_from_file(matrix_t A, char path[]) {
    std::ifstream file;
    file.open(path);

    for(size_t i = 0; i < A.rows * A.cols; i++)
        file >> A.ptr[i];

    file.close();
}

void matrix_to_file(matrix_t A, char path[]) {
    std::ofstream output(path);
    int n = A.rows;

    for (size_t i = 0; i < A.rows; ++i) {
        for (size_t j = 0; j < A.cols; ++j) {
            output << A.ptr[n*i +j] << " \t"; 
        }
        output << "\n";
    }
    output.close();
}
