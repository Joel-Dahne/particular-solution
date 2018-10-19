#ifndef SIGMA
#define SIGMA

#include "geom.h"
#include "mpfr.h"
#include "arb.h"

void sigma(mpfr_t res, geom_t geom, points_t points, slong N, arb_t nu,
           slong prec);

void minimize_sigma(arb_t nu, geom_t geom, points_t points, slong N,
                    arb_t nu_enclosure, arb_t tol, slong prec);

void coefs_sigma(arb_ptr coefs, geom_t geom, points_t points, slong N, arb_t nu,
                 slong prec);

#endif
