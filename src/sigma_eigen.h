#ifndef SIGMA_EIGEN
#define SIGMA_EIGEN

#include "mpfr.h"

void sigma_eigen(mpfr_t res, mpfr_t *A, int boundary, int rows, int N);

void coefs_sigma_eigen(mpfr_t *coefs, mpfr_t *A, int boundary, int rows,
                       int N);

#endif
