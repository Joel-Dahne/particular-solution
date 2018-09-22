#ifndef SIGMA
#define SIGMA

#include "mpfr.h"
#include "arb.h"

void sigma(mpfr_t res, points_t points, slong N, arb_t nu, arb_t mu0,
           int (*index)(int), slong prec);

void minimize_sigma(arb_t nu, points_t points, slong N, arb_t nu_enclosure,
                    arb_t mu0, arb_t tol, int (*index)(int), slong prec);

void coefs_sigma(arb_ptr coefs, points_t points, slong N, arb_t nu,
                 arb_t mu0, int (*index)(int), slong prec);

#endif
