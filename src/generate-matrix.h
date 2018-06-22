#ifndef GENERATE_MATRIX
#define GENERATE_MATRIX

#include "mpfr.h"

void scale_norm(mpfr_t scaling, mpfr_t theta_bound, mpfr_t nu, mpfr_t mu);

void generate_matrix(mpfr_t *A, mpfr_t *thetas, mpfr_t *phis, mpfr_t *scaling,
                     int len, int N, mpfr_t nu, mpfr_t mu0);

void eigenfunction(mpfr_t *res, mpfr_t *coefs, mpfr_t *thetas,
                   mpfr_t *phis, int len, int N, mpfr_t nu, mpfr_t mu0);

#endif
