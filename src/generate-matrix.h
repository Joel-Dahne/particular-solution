#ifndef GENERATE_MATRIX
#define GENERATE_MATRIX

#include "geom.h"
#include "arb.h"

void generate_matrix(mpfr_t *A, geom_t geom, points_t points, slong N, arb_t nu,
                     slong prec);

void eigenfunction(arb_ptr res, geom_t geom, arb_ptr coefs, points_t points,
                   slong N, arb_t nu, slong prec);

#endif
