#pragma once

/*
 * Коэффициенты cos sin вращения для ij элемента.
 */
struct svd_2x2_params {
    double d1, d2;
    double c1, s1;
    double c2, s2;
    double k;
};

/*
 *
 * @param w элемент в позиции (i,i).
 * @param x элемент в позиции (i,j).
 * @param y элемент в позиции (j,i).
 * @param z элемент в позиции (j,j).
 * @return коэффициенты матрицы 2x2.
 */
struct svd_2x2_params nsvd(double w, double x, double y, double z);
