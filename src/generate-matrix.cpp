#include "tools.h"

#include "arb.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_calc.h"

void generate_matrix(mpfr_t *A, struct Points points, int N, mpfr_t nu_mpfr,
                     mpfr_t mu0_mpfr, int (*index)(int)) {
  arb_t theta, phi, nu, mu0, mu, tmp, res;
  fmpz_t mu_int;
  slong prec, prec_local, n;

  arb_init(theta);
  arb_init(phi);
  arb_init(nu);
  arb_init(mu0);
  arb_init(mu);
  arb_init(tmp);
  arb_init(res);

  fmpz_init(mu_int);

  prec = mpfr_get_default_prec();

  arf_set_mpfr(arb_midref(nu), nu_mpfr);
  arf_set_mpfr(arb_midref(mu0), mu0_mpfr);

  n = points.boundary + points.interior;

  for (slong i = 0; i < n; i++) {
    arf_set_mpfr(arb_midref(theta), points.thetas[i]);
    arf_set_mpfr(arb_midref(phi), points.phis[i]);
    for (slong j = 0; j < N; j++) {
      arb_mul_si(mu, mu0, index(j), prec);
      if (arb_get_unique_fmpz(mu_int, mu)) {
        arb_set_fmpz(mu, mu_int);
      }
      arb_cos(tmp, theta, prec);
      prec_local = prec;
      do {
        arb_hypgeom_legendre_p(res, nu, mu, tmp, 0, prec_local);
        prec_local *=2;
      } while (!arb_is_finite(res));

      arb_mul(tmp, mu, phi, prec);
      arb_sin(tmp, tmp, prec);

      arb_mul(res, res, tmp, prec);

      arf_get_mpfr(A[j*n + i], arb_midref(res), MPFR_RNDN);
    }
  }

  arb_clear(theta);
  arb_clear(phi);
  arb_clear(nu);
  arb_clear(mu0);
  arb_clear(mu);
  arb_clear(tmp);
  arb_clear(res);

  fmpz_clear(mu_int);
}

void eigenfunction(mpfr_t *res, mpfr_t *coefs_mpfr, struct Points points, int N,
                   mpfr_t nu_mpfr, mpfr_t mu0_mpfr, int (*index)(int)) {
  arb_ptr coefs;
  arb_t theta, phi, nu, mu0, mu, tmp, term, sum;
  fmpz_t mu_int;
  slong prec, prec_local;

  coefs = _arb_vec_init(N);
  arb_init(theta);
  arb_init(phi);
  arb_init(nu);
  arb_init(mu0);
  arb_init(mu);
  arb_init(tmp);
  arb_init(term);
  arb_init(sum);

  fmpz_init(mu_int);

  prec = mpfr_get_default_prec();

  for (slong j = 0; j < N; j++) {
    arf_set_mpfr(arb_midref(coefs + j), coefs_mpfr[j]);
  }

  arf_set_mpfr(arb_midref(nu), nu_mpfr);
  arf_set_mpfr(arb_midref(mu0), mu0_mpfr);

  for (slong i = 0; i < points.boundary + points.interior; i++) {
    arf_set_mpfr(arb_midref(theta), points.thetas[i]);
    arf_set_mpfr(arb_midref(phi), points.phis[i]);

    arb_zero(sum);

    for (slong j = 0; j < N; j++) {
      arb_mul_si(mu, mu0, index(j), prec);
      if (arb_get_unique_fmpz(mu_int, mu)) {
        arb_set_fmpz(mu, mu_int);
      }

      arb_cos(tmp, theta, prec);
      prec_local = prec;
      do {
        arb_hypgeom_legendre_p(term, nu, mu, tmp, 0, prec_local);
        prec_local *=2;
      } while (!arb_is_finite(term));

      arb_mul(tmp, mu, phi, prec);
      arb_sin(tmp, tmp, prec);

      arb_mul(term, term, tmp, prec);
      arb_mul(term, term, coefs + j, prec);

      arb_add(sum, sum, term, prec);
    }

    arf_get_mpfr(res[i], arb_midref(sum), MPFR_RNDN);
  }

  _arb_vec_clear(coefs, N);
  arb_clear(theta);
  arb_clear(phi);
  arb_clear(nu);
  arb_clear(mu0);
  arb_clear(mu);
  arb_clear(tmp);
  arb_clear(term);
  arb_clear(sum);

  fmpz_clear(mu_int);
}
