#ifndef __MY_BLAS__
#define __MY_BLAS__

void xpay(const double, const double *, double *, const int);
void axpy(const double, const double *, double *, const int);
void dot(double *, const double *, const double *, const int);
void norm(double *, const double *, const int);

#endif
