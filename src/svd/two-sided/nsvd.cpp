#include "nsvd.hpp"
#include <math.h>
#include "../../utils/util.hpp"

struct svd_2x2_params nsvd(double w, double x, double y, double z) {
    bool flag = false;
    double c = 0.0, c1 = 0.0, c2 = 0.0;
    double s = 0.0, s1 = 0.0, s2 = 0.0;
    double k = 1.0;

    if (y == 0.0 && z == 0.0) {
        y = x;
        x = 0.0;
        flag = true;
    }

    double m1 = w + z;
    double m2 = x - y;

    if (!is_normalized(m2)) {
        c = 1.0;
        s = 0.0;
    } else {
        double r = m1 / m2;
        s = sign(r) / sqrt(1.0 + r * r);
        c = s * r;
    }

    m1 = s * (x + y) + c * (z - w);
    m2 = 2.0 * (c * x - s * z);

    if (!is_normalized(m2)) {
        c2 = 1.0;
        s2 = 0.0;
    } else {
        double r = m1 / m2;
        double t = sign(r) / (fabs(r) + sqrt(1.0 + r * r));
        c2 = 1.0 / sqrt(1.0 + t * t);
        s2 = c2 * t;
    }

    c1 = c2 * c - s2 * s;
    s1 = s2 * c + c2 * s;
    double d1 = c1 * (w * c2 - x * s2) - s1 * (y * c2 - z * s2);
    double d2 = s1 * (w * s2 + x * c2) + c1 * (y * s2 + z * c2);

    if (flag) {
        c2 = c1;
        s2 = s1;
        c1 = 1.0;
        s1 = 0.0;
    }
    if (fabs(d2) > fabs(d1)) {
        double r = c1;
        c1 = -s1;
        s1 = r;
        r = c2;
        c2 = -s2;
        s2 = r;
        r = d1;
        d1 = d2;
        d2 = r;
    }

    if (d1 < 0) {
        d1 = -d1;
        c1 = -c1;
        s1 = -s1;
        k = -k;
    }
    if (d2 < 0) {
        d2 = -d2;
        k = -k;
    }

    return {d1, d2, c1, s1, c2, s2, k};
}
