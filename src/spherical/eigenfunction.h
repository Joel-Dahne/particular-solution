#ifndef EIGENFUNCTION
#define EIGENFUNCTION

#include "geom.h"
#include "arb.h"

void eigenfunction_basis(arb_t res, arb_t theta, arb_t phi, arb_t nu, arb_t mu,
                         slong prec);

void eigenfunction_basis_series(arb_ptr res, arb_ptr z, arb_ptr phi, arb_t nu,
                                arb_ptr mu, slong n, slong prec);

void eigenfunction(arb_t res, geom_t geom, arb_ptr coefs, slong N, arb_t theta,
                   arb_t phi, arb_t nu, slong vertex, slong prec);

void eigenfunction_series(arb_ptr res, geom_t geom, arb_ptr coefs, slong N,
                          arb_ptr z, arb_ptr phi, arb_t nu, slong vertex,
                          slong n, slong prec);

#endif
