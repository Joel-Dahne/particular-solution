#ifndef NORM
#define NORM

#include "geom.h"
#include "arb.h"

void integral_norm(arb_t norm, geom_t geom, arb_ptr coefs, int N, arb_t mu,
                   slong vertex, slong prec);

#endif
