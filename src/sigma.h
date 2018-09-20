#ifndef SIGMA
#define SIGMA

#include "mpfr.h"
#include "arb.h"

void sigma(mpfr_t res, points_t points, int N, arb_t nu, arb_t mu0,
           int (*index)(int));

void minimize_sigma(arb_t nu, points_t points, int N, arb_t nu_enclosure,
                    arb_t mu0, arb_t tol, int (*index)(int));

void coefs_sigma(arb_ptr coefs, points_t points, int N, arb_t nu,
                 arb_t mu0, int (*index)(int));

#endif
