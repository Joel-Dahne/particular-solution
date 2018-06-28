#ifndef GENERATE_MATRIX
#define GENERATE_MATRIX

#include "mpfr.h"

void scale_norm(mpfr_t scaling, mpfr_t theta_bound, mpfr_t nu, mpfr_t mu);

void generate_matrix(mpfr_t *A, struct Points points, mpfr_t *scaling, int N,
                     mpfr_t nu, mpfr_t mu0, int (*index)(int));

void eigenfunction(mpfr_t *res, mpfr_t *coefs, struct Points points, int N,
                   mpfr_t nu, mpfr_t mu0, int (*index)(int));

#endif
