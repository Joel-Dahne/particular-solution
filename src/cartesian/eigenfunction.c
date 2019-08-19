#include "eigenfunction.h"

#include "arb_hypgeom.h"

/* Compute the value of the basis functions used in the expansion for
 * the eigenfunction at the point given by angle and radii with the
 * paremeters nu and mu.
 *
 * nu: square root of the eigenvalue
 * mu: corresponds to nu in the standard notation
 */
void
eigenfunction_basis(arb_t res, arb_t r, arb_t theta, arb_t lambda, arb_t nu,
                    slong prec)
{
  arb_t tmp;

  arb_init(tmp);

  arb_sqrt(tmp, lambda, prec);
  arb_mul(tmp, tmp, r, prec);
  arb_hypgeom_bessel_j(tmp, nu, tmp, prec);

  arb_mul(res, nu, theta, prec);
  arb_sin(res, res, prec);

  arb_mul(res, res, tmp, prec);

  arb_clear(tmp);
}

/* TODO: Implement this*/
void
eigenfunction_basis_series(arb_ptr res, arb_ptr z, arb_ptr phi, arb_t nu,
                           arb_ptr mu, slong n, slong prec)
{

}

void
eigenfunction(arb_t res, geom_t geom, arb_ptr coefs, slong N, arb_t r,
                   arb_t theta, arb_t lambda, slong vertex, slong prec)
{
  arb_t nu, term;

  arb_init(nu);
  arb_init(term);

  arb_zero(res);

  for (slong i = 0; i < N; i++)
  {
    geom_get_mu(nu, geom, vertex, i, prec);

    eigenfunction_basis(term, r, theta, lambda, nu, prec);

    arb_addmul(res, term, coefs + i, prec);
  }

  arb_clear(nu);
  arb_clear(term);
}

/* TODO: Implement this*/
void
eigenfunction_series(arb_ptr res, geom_t geom, arb_ptr coefs, slong N,
                     arb_ptr z, arb_ptr phi, arb_t nu, slong vertex,
                     slong n, slong prec)
{

}
