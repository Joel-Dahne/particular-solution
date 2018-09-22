#ifndef GENERATE_MATRIX
#define GENERATE_MATRIX

#include "geom.h"
#include "arb.h"

void generate_matrix(mpfr_t *A, points_t points, slong N, arb_t nu, arb_t mu0,
                     int (*index)(int), slong prec);

void eigenfunction(mpfr_t *res, arb_ptr coefs, points_t points, slong N,
                   arb_t nu, arb_t mu0, int (*index)(int), slong prec);

#endif
