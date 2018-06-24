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

void scale_norm(mpfr_t mpfr_scaling, mpfr_t mpfr_theta_bound,
                mpfr_t mpfr_nu, mpfr_t mpfr_mu) {
  acb_ptr params;
  acb_t scaling, zero, theta_bound;
  mag_t tol;
  slong prec;

  params = _acb_vec_init(3);
  acb_init(scaling);
  acb_init(zero);
  acb_init(theta_bound);
  mag_init(tol);

  prec = mpfr_get_default_prec();

  // The parameters are [nu, mu, tmp]
  arf_set_mpfr(arb_midref(acb_realref(params)), mpfr_nu);
  arf_set_mpfr(arb_midref(acb_realref(params + 1)), mpfr_mu);

  arb_set_str(acb_realref(zero), "1e-1", prec);
  arf_set_mpfr(arb_midref(acb_realref(theta_bound)), mpfr_theta_bound);

  mag_zero(tol);
  do {
    acb_calc_integrate(scaling, integrand, params, zero, theta_bound, prec, tol,
                       NULL, prec);
    prec *= 2;
  } while (!acb_is_finite(scaling));

  arf_get_mpfr(mpfr_scaling, arb_midref(acb_realref(scaling)), MPFR_RNDN);
  mpfr_rec_sqrt(mpfr_scaling, mpfr_scaling, MPFR_RNDN);

  _acb_vec_clear(params, 3);
  acb_clear(scaling);
  acb_clear(zero);
  acb_clear(theta_bound);
  mag_clear(tol);
}

void generate_matrix(mpfr_t *A_arr, mpfr_t *thetas, mpfr_t *phis,
                     mpfr_t *scaling, int len, int N, mpfr_t nu, mpfr_t mu0) {
  arb_t theta, phi, arb_nu, arb_mu0, tmp, tmp2, res;
  slong prec, prec_local;

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

  for (slong i = 0; i < len; i++) {
    arf_set_mpfr(arb_midref(theta), thetas[i]);
    arf_set_mpfr(arb_midref(phi), phis[i]);
    for (slong j = 0; j < N; j++) {
      arb_mul_si(tmp, arb_mu0, 2*j + 1, prec);

      arb_cos(tmp2, theta, prec);
      prec_local = prec;
      do {
        arb_hypgeom_legendre_p(res, arb_nu, tmp, tmp2, 0, prec_local);
        prec_local *=2;
      } while (!arb_is_finite(res));

      arb_mul(tmp, tmp, phi, prec);
      arb_sin(tmp, tmp, prec);

      arb_mul(res, res, tmp, prec);

      arf_get_mpfr(A_arr[j*len + i], arb_midref(res), MPFR_RNDN);
      mpfr_mul(A_arr[j*len + i], A_arr[j*len + i], scaling[j], MPFR_RNDN);
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

void eigenfunction(mpfr_t *res, mpfr_t *coefs, mpfr_t *thetas,
                   mpfr_t *phis, int len, int N, mpfr_t nu, mpfr_t mu0) {
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

  for (slong i = 0; i < len; i++) {
    arf_set_mpfr(arb_midref(theta), thetas[i]);
    arf_set_mpfr(arb_midref(phi), phis[i]);

    arb_zero(sum);

    for (slong j = 0; j < N; j++) {
      arb_mul_si(tmp, arb_mu0, 2*j + 1, prec);

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
