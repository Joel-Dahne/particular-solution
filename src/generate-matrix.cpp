#include "generate-matrix.h"

#include "arb_hypgeom.h"

static void
eval_eigenfunction(arb_t res, arb_t theta, arb_t phi, arb_t nu, arb_t mu,
                   slong prec)
{
  arb_t tmp;
  slong prec_local;

  arb_init(tmp);

  arb_cos(tmp, theta, prec);
  prec_local = prec;
  do {
    arb_hypgeom_legendre_p(res, nu, mu, tmp, 0, prec_local);
    prec_local *=2;
  } while (!arb_is_finite(res));

  arb_mul(tmp, mu, phi, prec);
  arb_sin(tmp, tmp, prec);

  arb_mul(res, res, tmp, prec);

  arb_clear(tmp);
}

void generate_matrix(mpfr_t *A, geom_t geom, points_t points, slong N, arb_t nu,
                     slong prec) {
  arb_t mu, res;
  slong n;

  arb_init(mu);
  arb_init(res);

  n = points->boundary + points->interior;

  for (slong i = 0; i < n; i++) {
    for (slong j = 0; j < N; j++) {
      geom_get_mu(mu, geom, 0, j, prec);

      eval_eigenfunction(res, points->thetas + i, points->phis + i, nu, mu,
                         prec);

      arf_get_mpfr(A[j*n + i], arb_midref(res), MPFR_RNDN);
    }
  }

  arb_clear(mu);
  arb_clear(res);
}

void eigenfunction(arb_ptr res, geom_t geom, arb_ptr coefs, points_t points,
                   slong N, arb_t nu, slong prec) {
  arb_t mu, term;

  arb_init(mu);
  arb_init(term);

  for (slong i = 0; i < points->boundary + points->interior; i++) {
    arb_zero(res + i);
    for (slong j = 0; j < N; j++) {
      geom_get_mu(mu, geom, 0, j, prec);

      eval_eigenfunction(term, points->thetas + i, points->phis + i, nu, mu, prec);

      arb_addmul(res + i, term, coefs + j, prec);
    }
  }

  arb_clear(mu);
  arb_clear(term);
}
