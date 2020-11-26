#include "svd_subprocedure.hpp"
#include <math.h>
#include "../../utils/matrix.hpp"
#include "nsvd.hpp"

size_t svd_subprocedure(struct matrix_t Bmat, struct matrix_t Umat, struct matrix_t Vmat) {
    size_t iter = 0;  //число повторений цикла развертки
    size_t n = Bmat.rows; //размер матрицы
    const double tol = 1e-15; //точность предела сходимости
    double* B = Bmat.ptr;
    double* U = Umat.ptr;
    double* V = Vmat.ptr;
    double norm = 0.0;      //норма Фробениуса матрицы B
    double off_norm = 0.0;  //норма Фробениуса только недиагональных элементов матрицы B

    //матрицы U и V инициализируются как единичные матрицы
    matrix_identity(Umat);
    matrix_identity(Vmat);
    //вычисление нормы Фробениуса, и нормы недиагональных элементов
    matrix_frobenius(Bmat, &norm, &off_norm);

    while (sqrt(off_norm) > tol * sqrt(norm)) {
        for (size_t i = 0; i < n - 1; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                double bii = B[n * i + i];
                double bij = B[n * i + j];
                double bji = B[n * j + i];
                double bjj = B[n * j + j];

                //вычислить коэффициенты cos sin для элемента bij
                //методом NSVD - https://maths-people.anu.edu.au/~brent/pd/rpb080i.pdf (стр. 12)
                struct svd_2x2_params cf = nsvd(bii, bij, bji, bjj);

                //выполнение операции J^T*B*J, т.е. обновление строк и столбцов
                //J - матрица плоских вращений (матрица Гивенса) с cos sin на местах i j
                //первый цикл обновляет строки i j
                for (size_t k = 0; k < n; k++) {
                    double b_ik = B[n * i + k];
                    double b_jk = B[n * j + k];

                    double left = cf.c1 * b_ik - cf.s1 * b_jk;
                    double right = cf.s1 * cf.k * b_ik + cf.c1 * cf.k * b_jk;

                    B[n * i + k] = left;
                    B[n * j + k] = right;
                }

                //второй цикл обновляет столбцы i j
                for (size_t k = 0; k < n; k++) {
                    double b_ki = B[n * k + i];
                    double b_kj = B[n * k + j];

                    double left = cf.c2 * b_ki - cf.s2 * b_kj;
                    double right = cf.s2 * b_ki + cf.c2 * b_kj;

                    B[n * k + i] = left;
                    B[n * k + j] = right;
                }

                //последние два цикла обновляют столбцы i j у матриц U V 
                for (size_t k = 0; k < n; k++) {
                    double u_ki = U[n * k + i];
                    double u_kj = U[n * k + j];

                    double left = cf.c1 * u_ki - cf.s1 * u_kj;
                    double right = cf.s1 * cf.k * u_ki + cf.c1 * cf.k * u_kj;

                    U[n * k + i] = left;
                    U[n * k + j] = right;
                }

                for (size_t k = 0; k < n; k++) {
                    double v_ki = V[n * k + i];
                    double v_kj = V[n * k + j];

                    double left = cf.c2 * v_ki - cf.s2 * v_kj;
                    double right = cf.s2 * v_ki + cf.c2 * v_kj;

                    V[n * k + i] = left;
                    V[n * k + j] = right;
                }
            }
        }

        matrix_frobenius(Bmat, &norm, &off_norm);
        iter += 1;
    }

    return iter;
}
