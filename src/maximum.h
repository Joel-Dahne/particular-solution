#ifndef MAXIMUM
#define MAXIMUM

#include "geom.h"
#include "arb.h"

void maximize_series(arb_t max, geom_t geom, arb_t t, arb_ptr coefs, slong N,
                     arb_t nu, slong n, slong vertex, slong prec);

void maximize(arb_t max, geom_t geom, arb_ptr coefs, slong N, arb_t nu,
              slong vertex, slong prec);

#endif
