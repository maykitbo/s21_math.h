#ifndef S21_MATH_H_
#define S21_MATH_H_

#define s21_M_PI 3.14159265358979323846
#define s21_M_E 2.7182818284590452353602874713527
#define s21_ln_0_5 -0.6931471805599453
#define s21_ln_2 0.6931471805599453
#define s21_INFINITY 1/0.0
#define s21_NAN 0.0/0.0
#define s21_sqrt_10 3.1622776601683795
#define s21_tan_05pi 16331239353195370.000000

int s21_abs(int x);
long double s21_fabs(double x);
long double s21_ceil(double x);
long double s21_floor(double x);
long double s21_pow(double x, double n);
long double s21_sin(double x);
long double s21_cos(double x);
long double s21_tan(double x);
long double s21_log(double x);
long double s21_asin(double x);
long double s21_acos(double x);
long double s21_atan(double x);
long double s21_exp(double x);
long double s21_sqrt(double x);
long double s21_fmod(double x, double y);
int s21_isnan(double x);
int s21_isinf(double x);

#endif