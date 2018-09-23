#include "generate-matrix.h"

#include "arb_hypgeom.h"

void generate_matrix(mpfr_t *A, points_t points, slong N, arb_t nu,
                     arb_t mu0, int (*index)(int), slong prec) {
  arb_t mu, tmp, res;
  fmpz_t mu_int;
  slong prec_local, n;

  arb_init(mu);
  arb_init(tmp);
  arb_init(res);

  fmpz_init(mu_int);

  n = points->boundary + points->interior;

  for (slong i = 0; i < n; i++) {
    for (slong j = 0; j < N; j++) {
      arb_mul_si(mu, mu0, index(j), prec);
      if (arb_get_unique_fmpz(mu_int, mu)) {
        arb_set_fmpz(mu, mu_int);
      }
      arb_cos(tmp, points->thetas + i, prec);
      prec_local = prec;
      do {
        arb_hypgeom_legendre_p(res, nu, mu, tmp, 0, prec_local);
        prec_local *=2;
      } while (!arb_is_finite(res));

      arb_mul(tmp, mu, points->phis + i, prec);
      arb_sin(tmp, tmp, prec);

      arb_mul(res, res, tmp, prec);

      arf_get_mpfr(A[j*n + i], arb_midref(res), MPFR_RNDN);
    }
  }

  arb_clear(mu);
  arb_clear(tmp);
  arb_clear(res);

  fmpz_clear(mu_int);
}

void eigenfunction(arb_ptr res, arb_ptr coefs, points_t points, slong N,
                   arb_t nu, arb_t mu0, int (*index)(int), slong prec) {
  arb_t mu, tmp, term;
  fmpz_t mu_int;
  slong prec_local;

  arb_init(mu);
  arb_init(tmp);
  arb_init(term);

  fmpz_init(mu_int);

  for (slong i = 0; i < points->boundary + points->interior; i++) {
    arb_zero(res + i);
    for (slong j = 0; j < N; j++) {
      arb_mul_si(mu, mu0, index(j), prec);
      if (arb_get_unique_fmpz(mu_int, mu)) {
        arb_set_fmpz(mu, mu_int);
      }

      arb_cos(tmp, points->thetas + i, prec);
      prec_local = prec;
      do {
        arb_hypgeom_legendre_p(term, nu, mu, tmp, 0, prec_local);
        prec_local *=2;
      } while (!arb_is_finite(term));

      arb_mul(tmp, mu, points->phis + i, prec);
      arb_sin(tmp, tmp, prec);

      arb_mul(term, term, tmp, prec);

      arb_addmul(res + i, term, coefs + j, prec);
    }
  }

  arb_clear(mu);
  arb_clear(tmp);
  arb_clear(term);

  fmpz_clear(mu_int);
}
