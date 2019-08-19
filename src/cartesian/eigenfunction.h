#ifndef EIGENFUNCTION
#define EIGENFUNCTION

#include "geom.h"
#include "arb.h"

void eigenfunction_basis(arb_t res, arb_t r, arb_t theta, arb_t lambda, arb_t nu,
                         slong prec);

void eigenfunction_basis_series(arb_ptr res, arb_t r, arb_t theta, arb_t lambda,
                                arb_t nu, slong n, slong prec);

void eigenfunction(arb_t res, geom_t geom, arb_ptr coefs, slong N, arb_t r,
                   arb_t theta, arb_t lambda, slong vertex, slong prec);

void eigenfunction_series(arb_ptr res, geom_t geom, arb_ptr coefs, slong N,
                          arb_ptr r, arb_ptr theta, arb_t lambda, slong vertex,
                          slong n, slong prec);

#endif
