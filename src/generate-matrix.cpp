#include "tools.h"

#include "arb.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_calc.h"

int integrand(acb_ptr out, const acb_t inp, void *param, slong order,
              slong prec) {
  acb_ptr pr;
  slong prec_local;

  pr = (acb_ptr)(param);

  if (order > 1)
    flint_abort();

  // We check intersection with the branch cut by checking if contains
  // 1, this is perhaps not optimal.
  if (order == 1 && arb_contains_si(acb_realref(inp), 1)) {
    acb_indeterminate(out);
    return 0;
  }

  acb_cos(out, inp, prec);

  prec_local = prec;
  do {
    acb_hypgeom_legendre_p(pr + 2, pr, pr + 1, out, 0, prec_local);
    prec_local *= 2;
  } while (!acb_is_finite(pr + 2) && prec_local <= 16*prec);

  acb_sqr(out, pr + 2, prec);

  acb_sin(pr + 2, inp, prec);

  acb_mul(out, out, pr + 2, prec);

  return 0;
}

void generate_matrix(mpfr_t *A_arr, struct Points points, int N, mpfr_t nu,
                     mpfr_t mu0, int (*index)(int)) {
  arb_t theta, phi, arb_nu, arb_mu0, tmp, tmp2, res;
  slong prec, prec_local, n;

  arb_init(theta);
  arb_init(phi);
  arb_init(arb_nu);
  arb_init(arb_mu0);
  arb_init(tmp);
  arb_init(tmp2);
  arb_init(res);

  prec = mpfr_get_default_prec();

  arf_set_mpfr(arb_midref(arb_nu), nu);
  arf_set_mpfr(arb_midref(arb_mu0), mu0);

  n = points.boundary + points.interior;

  for (slong i = 0; i < n; i++) {
    arf_set_mpfr(arb_midref(theta), points.thetas[i]);
    arf_set_mpfr(arb_midref(phi), points.phis[i]);
    for (slong j = 0; j < N; j++) {
      arb_mul_si(tmp, arb_mu0, index(j), prec);

      arb_cos(tmp2, theta, prec);
      prec_local = prec;
      do {
        arb_hypgeom_legendre_p(res, arb_nu, tmp, tmp2, 0, prec_local);
        prec_local *=2;
      } while (!arb_is_finite(res));

      arb_mul(tmp, tmp, phi, prec);
      arb_sin(tmp, tmp, prec);

      arb_mul(res, res, tmp, prec);

      arf_get_mpfr(A_arr[j*n + i], arb_midref(res), MPFR_RNDN);
    }
  }

  arb_clear(theta);
  arb_clear(phi);
  arb_clear(arb_nu);
  arb_clear(arb_mu0);
  arb_clear(tmp);
  arb_clear(tmp2);
  arb_clear(res);
}

void eigenfunction(mpfr_t *res, mpfr_t *coefs, struct Points points, int N,
                   mpfr_t nu, mpfr_t mu0, int (*index)(int)) {
  arb_ptr arb_coefs;
  arb_t theta, phi, arb_nu, arb_mu0, tmp, tmp2, term, sum;
  slong prec, prec_local;

  arb_coefs = _arb_vec_init(N);
  arb_init(theta);
  arb_init(phi);
  arb_init(arb_nu);
  arb_init(arb_mu0);
  arb_init(tmp);
  arb_init(tmp2);
  arb_init(term);
  arb_init(sum);

  prec = mpfr_get_default_prec();

  for (slong j = 0; j < N; j++) {
    arf_set_mpfr(arb_midref(arb_coefs + j), coefs[j]);
  }

  arf_set_mpfr(arb_midref(arb_nu), nu);
  arf_set_mpfr(arb_midref(arb_mu0), mu0);

  for (slong i = 0; i < points.boundary + points.interior; i++) {
    arf_set_mpfr(arb_midref(theta), points.thetas[i]);
    arf_set_mpfr(arb_midref(phi), points.phis[i]);

    arb_zero(sum);

    for (slong j = 0; j < N; j++) {
      arb_mul_si(tmp, arb_mu0, index(j), prec);

      arb_cos(tmp2, theta, prec);
      prec_local = prec;
      do {
        arb_hypgeom_legendre_p(term, arb_nu, tmp, tmp2, 0, prec_local);
        prec_local *=2;
      } while (!arb_is_finite(term));

      arb_mul(tmp, tmp, phi, prec);
      arb_sin(tmp, tmp, prec);

      arb_mul(term, term, tmp, prec);
      arb_mul(term, term, arb_coefs + j, prec);

      arb_add(sum, sum, term, prec);
    }

    arf_get_mpfr(res[i], arb_midref(sum), MPFR_RNDN);
  }

  _arb_vec_clear(arb_coefs, N);
  arb_clear(theta);
  arb_clear(phi);
  arb_clear(arb_nu);
  arb_clear(arb_mu0);
  arb_clear(tmp);
  arb_clear(tmp2);
  arb_clear(term);
  arb_clear(sum);
}
