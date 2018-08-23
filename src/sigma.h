#ifndef SIGMA
#define SIGMA

#include "mpfr.h"

void sigma(mpfr_t res, points_t points, int N, mpfr_t nu, mpfr_t mu0,
           int (*index)(int));

void minimize_sigma(mpfr_t nu, points_t points, int N, mpfr_t nu_low,
                    mpfr_t nu_upp, mpfr_t mu0, mpfr_t tol, int (*index)(int));

void coefs_sigma(mpfr_t *coefs_mpfr, points_t points, int N, mpfr_t nu,
                 mpfr_t mu0, int (*index)(int));

#endif
