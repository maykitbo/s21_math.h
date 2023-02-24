#include "s21_math.h"

static long double s21_floor_ceil(double x, int c, int f);
static long double reduction_sin_cos(double x);
static int check_input_pow(double *x, double n, int *znak_x, long double *result);
static long double n_not_integer(double x, double n, int *znak_n, int *znak_x, long double result);

long double s21_pow(double x, double n) {
    long double result = 1;
    int znak_x = 1;
    if (check_input_pow(&x, n, &znak_x, &result)) {
        int znak_n = 1;
        if (s21_ceil(n) == n) {
            if (n < 0) {
                n *= -1;
                znak_n = -1;
            }
            for (int k = 0; k < n; k++) result *= x;
        } else {
            result = n_not_integer(x, n, &znak_n, &znak_x, result);
        }
        if (znak_n == -1) result = 1 / result;
    }
    return result * znak_x;
}

static long double n_not_integer(double x, double n, int *znak_n, int *znak_x, long double result) {
    int n_big = 1;
    long double a = 1;
    if (x < 0) {
        x *= -1;
        if (s21_floor(n / 2) != n / 2) *znak_x = -1;
    }
    if (x > 2) {
        long double pow_2_n = s21_pow(2, n);
        while (x > 2) {
            result *= pow_2_n;
            x /= 2;
        }
        if (!s21_isinf(result)) result *= s21_pow(x, n);
    } else if (x < 1) {
        long double pow_2_n = s21_pow(2, n);
        while (x < 1) {
            result /= pow_2_n;
            x *= 2;
        }
        if (!s21_isinf(result)) result *= s21_pow(x, n);
    } else {
        a = 1;
        x -= 1;
        if (n < 0) {
            n *= -1;
            *znak_n = -1;
        }
        if (n < 1 && n > 0) {
            n += 1;
            n_big = 0;
        }
        long double mnoj = (n * x) * result / a;
        result += mnoj;
        long double buff = 0;
        for (int k = 2; k < 1e7; k++) {
            mnoj *= (x * (n - k + 1)) / (k * a);
            result += mnoj;
            if (result == buff) break;
            buff = result;
        }
    }
    if (n_big == 0) result /= x + a;
    return result;
}

static int check_input_pow(double *x, double n, int *znak_x, long double *result) {
    int res = 0;
    if (n == 0 || *x == 1) {
        *result = 1;
    } else if (s21_isnan(n)) {
        *result = s21_NAN;
    } else if (*x == -1) {
        if (s21_ceil(n / 2) == n / 2) *result = 1;
        else
            *result = -1;
    } else if (*x == 0) {
        if (n < 0) *result = s21_INFINITY;
        else
            *result = 0;
    } else if (s21_isnan(*x)) {
        *result = s21_NAN;
    } else if (*x == s21_INFINITY) {
        if (n > 0) *result = s21_INFINITY;
        else
            *result = 0;
    } else if (*x == -s21_INFINITY) {
        if (n < 0) *result = 0;
        else if (s21_ceil(n / 2) == n / 2) *result = s21_INFINITY;
        else
            *result = -s21_INFINITY;
    } else if (n == s21_INFINITY) {
        if (*x > 1) *result = n;
        else if (*x < 1 && *x > -1) *result = 0;
        else if (*x < -1 && s21_ceil(n / 2) != n / 2) *result = -s21_INFINITY;
        else
            *result = s21_INFINITY;
    } else if (n == -s21_INFINITY) {
        if (*x > 1 || *x < -1) *result = 0;
        else if (*x < 0) *result = -s21_INFINITY;
        else
            *result = s21_INFINITY;
    } else if (n == 1) {
        if (*x < 0) {
            *znak_x = -1;
            *x = -1 * (*x);
        }
        *result = *x;
    } else if (n == -1) {
        if (*x < 0) {
            *znak_x = -1;
            *x = -1 * (*x);
        }
        *result = 1 / (*x);
    } else if (*x < 0 && s21_floor(n) != n) {
        *result = -s21_NAN;
    } else {
        res = 1;
    }
    return res;
}

int s21_isnan(double x) {
    int result = 0;
    if (!(x >= 0) && !(x < 0)) result = 1;
    return result;
}

int s21_isinf(double x) {
    int result = 0;
    if (x == s21_INFINITY || x == -s21_INFINITY) result = 1;
    return result;
}

long double s21_sqrt(double x) {
    long double result = 1;
    if (x == - s21_INFINITY) result = s21_NAN;
    if (x != s21_INFINITY) {
        while (x > 10) {
            x /= 10;
            result *= s21_sqrt_10;
        }
    }
    result *= s21_pow(x, 0.5);
    return result;
}

int s21_abs(int x) {
    if (x < 0) x *= -1;
    return x;
}

long double s21_fabs(double x) {
    if (x < 0) x *= -1;
    return x;
}

long double s21_ceil(double x) {
    return s21_floor_ceil(x, 1, 0);
}

long double s21_floor(double x) {
    return s21_floor_ceil(x, 0, 1);
}

static long double s21_floor_ceil(double x, int c, int f) {
    double r = x;
    if (x != s21_INFINITY && x != -s21_INFINITY) {
        if (x > 1e6) {
            long double e = 1e308;
            for (int k = 308; k > 5; k--) {
                while (x > e) {
                    x -= e;
                }
                e /= 10;
            }
        }
        if ((int)x != x) {
            if (x > 0) r += (int)x - x + c;
            else
                r += (int)x - x - f;
        }
    }
    return r;
}

long double s21_fmod(double x, double y) {
    double result = s21_NAN;
    if (s21_isinf(x)) {
    } else if (s21_isinf(y)) {
        result = x;
    } else {
        int znak = 1;
        if (x < 0) {
            znak = -1;
            x *= -1;
        }
        y = s21_fabs(y);
        if (y != 0) result = (x - s21_floor(x / y) * y) * znak;
    }
    return result;
}

long double s21_sin(double x) {
    x = reduction_sin_cos(x);
    long double result = x, mnoj = x;
    if (!s21_isnan(x)) {
        long int k = 3;
        double xx = x * x;
        while (s21_fabs(mnoj) > 1e-7) {
            mnoj *= (-1 * xx) / (k * (k - 1));
            result += mnoj;
            k += 2;
        }
    }
    return result;
}

long double s21_cos(double x) {
    x = reduction_sin_cos(x);
    long double result = x;
    if (!s21_isnan(x)) {
        result = 1;
        long double mnoj = 1;
        long int k = 2;
        while (s21_fabs(mnoj) > 1e-7) {
            mnoj *= -1 * x * x / (k * (k - 1));
            result +=mnoj;
            k += 2;
        }
    }
    return result;
}

static long double reduction_sin_cos(double x) {
    if (s21_isinf(x)) x = -s21_NAN;
    if (x < - 2 * s21_M_PI || x > 2 * s21_M_PI)
        x -= s21_floor(x / (2 * s21_M_PI)) * 2 * s21_M_PI;
    return x;
}

long double s21_tan(double x) {
    long double result;
    if (s21_isinf(x)) result = s21_NAN;
    else if ((x + s21_M_PI / 2) / 2 * s21_M_PI == s21_ceil((x + s21_M_PI / 2) / 2 * s21_M_PI))
        result = -s21_tan_05pi;
    else if ((x + s21_M_PI / 2) / s21_M_PI == s21_ceil((x + s21_M_PI / 2) / s21_M_PI))
        result = s21_tan_05pi;
    else
        result = s21_sin(x) / s21_cos(x);
    return result;
}

long double s21_asin(double x) {
    long double result;
    if (x > -1 && x < 1) {
        result = x + x * x * x / 6;
        double l = x * x / 2;
        long double mnoj = 0.5 * x * x * x;
        long int k = 2;
        while (s21_fabs(mnoj) / (2 * k + 1) > 1e-12) {
            mnoj *= (2 * k - 1) * l / k;
            result += mnoj / (2 * k + 1);
            k++;
        }
    } else if (x == 1) {
        result = s21_M_PI / 2;
    } else if (x == -1) {
        result = - s21_M_PI / 2;
    } else {
        result = s21_NAN;
    }
    return result;
}

long double s21_acos(double x) {
    long double result;
    if (x >= -1 && x <= 1) result = s21_M_PI / 2 - s21_asin(x);
    else
        result = s21_NAN;
    return result;
}

long double s21_atan(double x) {
    long double result;
    if (x == s21_INFINITY) result = s21_M_PI / 2;
    else if (x == -s21_INFINITY) result = -s21_M_PI / 2;
    else
        result = s21_asin(x / s21_pow(x * x + 1, 0.5));
    return result;
}

long double s21_log(double x) {
    long double result;
    if (x < 0) {
        result = -s21_NAN;
    } else if (s21_isnan(x) || s21_isinf(x)) {
        result = x;
    } else if (x == s21_M_E) {
        result = 1;
    } else if (x == 0) {
        result = -s21_INFINITY;
    } else if (x == 2) {
        result = s21_ln_2;
    } else {
        int exp = 0;
        while (x > 2) {
            x /= 2;
            exp++;
        }
        while (x < 1) {
            x *= 2;
            exp--;
        }
        x--;
        result = x;
        double mnoj = x;
        long int k = 2;
        while (s21_fabs(mnoj / k) > 2e-7) {
            mnoj *= (-1) * (x);
            result += mnoj / k;
            k++;
        }
        result += exp * s21_ln_2;
    }
    return result;
}

long double s21_exp(double x) {
    long double result = 1;
    if (x == s21_INFINITY) {
        result = x;
    } else if (x == -s21_INFINITY) {
        result = 0;
    } else {
        long double mnoj = 1, znak_x = 0;
        if (x < 0) {
            x *= -1;
            znak_x = 1;
        }
        double r = s21_floor(x);
        x = x - s21_floor(x);
        long int k = 1;
        while (mnoj > 1e-15) {
            mnoj *= x;
            mnoj /= k;
            result += mnoj;
            k++;
        }
        for (double i = 0; i < r; i++) result *= s21_M_E;
        if (result > 1e308) result = s21_INFINITY;
        if (znak_x) result = 1 / result;
    }
    return result;
}
