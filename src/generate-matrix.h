#ifndef GENERATE_MATRIX
#define GENERATE_MATRIX

#include "geom.h"
#include "arb.h"

void generate_matrix(mpfr_t *A, points_t points, int N, mpfr_t nu, arb_t mu0,
                     int (*index)(int));

void eigenfunction(mpfr_t *res, mpfr_t *coefs, points_t points, int N,
                   mpfr_t nu, arb_t mu0, int (*index)(int));

#endif
