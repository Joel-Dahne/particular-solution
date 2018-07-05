#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_calc.h"

#include <iostream>

using namespace std;

static void angles_to_vectors_arb(arb_ptr v1, arb_ptr v2, arb_t theta_bound,
                                  arb_t critical_point, int angles_coefs[],
                                  slong prec) {
  arb_ptr angles, w;
  arb_t S, theta1, theta2, tmp, tmp2, tmp_dot;

  angles = _arb_vec_init(3);
  w = _arb_vec_init(3);
  arb_init(S);
  arb_init(theta1);
  arb_init(theta2);
  arb_init(tmp);
  arb_init(tmp2);
  arb_init(tmp_dot);

  for (slong i = 0; i < 3; i++) {
    arb_set_si(angles + i, angles_coefs[2*i]);
    arb_div_si(angles + i, angles + i, angles_coefs[2*i + 1], prec);
  }

  // Compute the sum of the angles
  arb_add(S, angles, angles + 1, prec);
  arb_add(S, S, angles + 2, prec);
  arb_div_si(S, S, 2, prec);

  /* Compute theta1 for v1*/
  arb_cos_pi(theta1, S, prec);

  arb_sub(tmp, S, angles + 1, prec);
  arb_cos_pi(tmp, tmp, prec);
  arb_mul(theta1, theta1, tmp, prec);

  arb_sin_pi(tmp, angles + 0, prec);
  arb_div(theta1, theta1, tmp, prec);

  arb_sin_pi(tmp, angles + 2, prec);
  arb_div(theta1, theta1, tmp, prec);

  arb_neg(theta1, theta1);
  arb_sqrt(theta1, theta1, prec);
  arb_asin(theta1, theta1, prec);
  arb_mul_si(theta1, theta1, 2, prec);

  /* Compute v1 from knowing theta1 */
  arb_sin(v1 + 0, theta1, prec);
  arb_set_si(v1 + 1, 0);
  arb_cos(v1 + 2, theta1, prec);

  /* Compute theta2 for v2*/
  arb_cos_pi(theta2, S, prec);

  arb_sub(tmp, S, angles + 2, prec);
  arb_cos_pi(tmp, tmp, prec);
  arb_mul(theta2, theta2, tmp, prec);

  arb_sin_pi(tmp, angles + 0, prec);
  arb_div(theta2, theta2, tmp, prec);

  arb_sin_pi(tmp, angles + 1, prec);
  arb_div(theta2, theta2, tmp, prec);

  arb_neg(theta2, theta2);
  arb_sqrt(theta2, theta2, prec);
  arb_asin(theta2, theta2, prec);
  arb_mul_si(theta2, theta2, 2, prec);

  /* Compute v2 from knowing theta2 */
  arb_sin(tmp, theta2, prec);
  arb_cos_pi(v2 + 0, angles + 0, prec);
  arb_mul(v2 + 0, v2 + 0, tmp, prec);
  arb_sin_pi(v2 + 1, angles + 0, prec);
  arb_mul(v2 + 1, v2 + 1, tmp, prec);
  arb_cos(v2 + 2, theta2, prec);

  /* Compute the critical point */
  _arb_vec_sub(w, v2, v1, 3, prec);
  /* dot(v, v) */
  _arb_vec_dot(tmp_dot, v1, v1, 3, prec);
  /* w_3*dot(v, v) */
  arb_mul(critical_point, w + 2, tmp_dot, prec);
  /* dot(v, w) */
  _arb_vec_dot(tmp_dot, v1, w, 3, prec);
  /* v_3*dot(v, w) */
  arb_mul(tmp, v1 + 2, tmp_dot, prec);
  /* w_3*dot(v, v) - v_3*dot(v, w)*/
  arb_sub(critical_point, tmp, critical_point, prec);
  /* w_3*dot(v, w) */
  arb_mul(tmp, w + 2, tmp_dot, prec);
  /* dot(w, w) */
  _arb_vec_dot(tmp_dot, w, w, 3, prec);
  /* v_3*dot(w, w) */
  arb_mul(tmp2, v1 + 2, tmp_dot, prec);
  /* w_3*dot(v, w) - v_3*dot(w, w)*/
  arb_sub(tmp, tmp, tmp2, prec);
  /* The critical point */
  arb_div(critical_point, critical_point, tmp, prec);

  /* If the critical point is in the interval [0, 1] compute the angle
   * corresponding to the critical point */
  arb_sub_si(tmp, critical_point, 1, prec);
  if (!(arb_contains_positive(critical_point) && arb_contains_negative(tmp)))
    arb_set_si(critical_point, 0);

  _arb_vec_scalar_mul(w, w, 3, critical_point, prec);
  _arb_vec_add(w, v1, w, 3, prec);
  _arb_vec_dot(tmp_dot, w, w, 3, prec);
  arb_sqrt(tmp_dot, tmp_dot, prec);
  /* Compute z for this point */
  arb_div(tmp, w + 2, tmp_dot, prec);
  arb_acos(tmp, tmp, prec);


  /* Compute the minimum angle */
  arb_min(theta_bound, theta1, theta2, prec);
  arb_min(theta_bound, theta_bound, tmp, prec);

  _arb_vec_clear(angles, 3);
  _arb_vec_clear(w, 3);
  arb_clear(S);
  arb_clear(theta1);
  arb_clear(theta2);
  arb_clear(tmp);
  arb_clear(tmp2);
  arb_clear(tmp_dot);
}

static int integrand(acb_ptr out, const acb_t inp, void *param_void, slong order,
                     slong prec) {
  acb_ptr param;
  slong prec_local;

  /* param = [nu, mu, tmp] */
  param = (acb_ptr)(param_void);

  if (order > 1)
    flint_abort();

  acb_cos(out, inp, prec);

  if (order == 1) {
    acb_sub_si(param + 2, out, 1, prec);
    if (!arb_is_nonpositive(acb_realref(param + 2))) {
      acb_indeterminate(out);
      return 0;
    }
  }

  prec_local = prec;
  do {
    acb_hypgeom_legendre_p(param + 2, param, param + 1, out, 0, prec_local);
    prec_local *= 2;
  } while (!acb_is_finite(param + 2) && prec_local <= 16*prec);

  acb_sqr(out, param + 2, prec);

  acb_sin(param + 2, inp, prec);

  acb_mul(out, out, param + 2, prec);

  return 0;
}

static void integral_norm(arb_t norm, arb_ptr coefs, int N, arb_t theta_bound,
                          arb_t nu, arb_t mu0, int (*index)(int), slong prec) {
  acb_t zero, theta_bound_c, tmp_c;
  acb_ptr params;
  arb_t integral_phi, mu, norm_part, tmp;
  mag_t tol;
  slong prec_local;

  params = _acb_vec_init(3);

  acb_init(zero);
  acb_init(theta_bound_c);
  acb_init(tmp_c);

  arb_init(integral_phi);
  arb_init(mu);
  arb_init(norm_part);
  arb_init(tmp);

  mag_init(tol);

  arb_set(acb_realref(params + 0), nu); /* params[0] = nu */

  arb_set_str(acb_realref(zero), "1e-1", prec);
  arb_set(acb_realref(theta_bound_c), theta_bound);

  /* Compute the integral corresponding to phi */
  arb_const_pi(integral_phi, prec);
  arb_div_si(integral_phi, integral_phi, 2, prec);
  arb_div(integral_phi, integral_phi, mu0, prec);
  arb_neg(integral_phi, integral_phi);

  for (slong i = 0; i < N; i++) {
    /* params[1] = mu */
    arb_mul_si(acb_realref(params + 1), mu0, index(i), prec);

    /* Compute the integral corresponding to theta */
    prec_local = prec;
    do {
      acb_calc_integrate(tmp_c, integrand, params, zero, theta_bound_c, prec, tol,
                         NULL, prec_local);
      prec_local *= 2;
    } while (!acb_is_finite(tmp_c));

    /* Multiply with the integral corresponding to phi */
    arb_mul(norm_part, acb_realref(tmp_c), integral_phi, prec);

    /* Multiply with the coefficient */
    arb_sqr(tmp, coefs + i, prec);
    arb_mul(norm_part, norm_part, tmp, prec);

    /* Add to the norm */
    arb_add(norm, norm, norm_part, prec);

  }
  arb_sqrt(norm, norm, prec);

  acb_clear(zero);
  acb_clear(theta_bound_c);
  acb_clear(tmp_c);

  _acb_vec_clear(params, 3);

  arb_clear(integral_phi);
  arb_clear(mu);
  arb_clear(norm_part);
  arb_clear(tmp);
}

/* It is important that the parametrization we have is good and
 * does not create larger balls then necessary. */
static void parametrization(arb_t z, arb_t phi, arf_t t_low, arf_t t_upp,
                             arb_ptr v1, arb_ptr v2, arb_t critical_point,
                             slong prec) {
  arb_ptr w;
  arb_t t, norm_w, z_low, z_upp, z_critical, tmp;
  arf_t phi_low, phi_upp, tmp2;
  int contains_critical;

  w = _arb_vec_init(3);
  arb_init(t);
  arb_init(norm_w);
  arb_init(z_low);
  arb_init(z_upp);
  arb_init(z_critical);
  arb_init(tmp);

  arf_init(phi_low);
  arf_init(phi_upp);
  arf_init(tmp2);

  /* Compute with the lower point */
  arf_set(arb_midref(t), t_low);

  _arb_vec_sub(w, v2, v1, 3, prec);
  _arb_vec_scalar_mul(w, w, 3, t, prec);
  _arb_vec_add(w, v1, w, 3, prec);

  _arb_vec_dot(norm_w, w, w, 3, prec);
  arb_sqrt(norm_w, norm_w, prec);
  _arb_vec_scalar_div(w, w, 3, norm_w, prec);

  arb_atan2(tmp, w + 1, w + 0, prec);
  arb_get_lbound_arf(phi_low, tmp, prec); // phi for the lower point
  arb_set(z_low, w + 2); // z for the lower point

  /* Compute with the upper point */
  arf_set(arb_midref(t), t_upp);

  _arb_vec_sub(w, v2, v1, 3, prec);
  _arb_vec_scalar_mul(w, w, 3, t, prec);
  _arb_vec_add(w, v1, w, 3, prec);

  _arb_vec_dot(norm_w, w, w, 3, prec);
  arb_sqrt(norm_w, norm_w, prec);
  _arb_vec_scalar_div(w, w, 3, norm_w, prec);

  arb_atan2(tmp, w + 1, w + 0, prec);
  arb_get_lbound_arf(phi_upp, tmp, prec); // phi for the upper point
  arb_set(z_upp, w + 2); // z for the upper point

  /* Compute with the critical point if relevant */
  arf_set(arb_midref(tmp), t_low);
  contains_critical = !arb_gt(tmp, critical_point); // t_low < critical_point
  arf_set(arb_midref(tmp), t_upp);
  contains_critical &= !arb_lt(tmp, critical_point); // t_low > critical_point

  if (contains_critical) {
    _arb_vec_sub(w, v2, v1, 3, prec);
    _arb_vec_scalar_mul(w, w, 3, t, prec);
    _arb_vec_add(w, v1, w, 3, prec);

    _arb_vec_dot(norm_w, w, w, 3, prec);
    arb_sqrt(norm_w, norm_w, prec);
    _arb_vec_scalar_div(w, w, 3, norm_w, prec);

    arb_atan2(tmp, w + 1, w + 0, prec);
    arb_set(z_critical, w + 2); // z for the critical point
  }

  /* phi is strictly increasing so is contained in the ball having
   * endpoints phi_low and phi_upp */
  arb_set_interval_arf(phi, phi_low, phi_upp, prec);

  /* z is not strictly increasing but has only one critical point. If
   * the critical point is not in the interval it is contained in the
   * ball having endpoint z_low and z_upp (not necessarily in that
   * order). If the critical point is contained in the interval we
   * take the ball containing all three values. */
  arb_union(z, z_low, z_upp, prec);
  if (contains_critical)
    arb_union(z, z, z_critical, prec);

  _arb_vec_clear(w, 3);
  arb_clear(t);
  arb_clear(norm_w);
  arb_clear(z_low);
  arb_clear(z_upp);
  arb_clear(z_critical);
  arb_clear(tmp);

  arf_clear(phi_low);
  arf_clear(phi_upp);
  arf_clear(tmp2);
}

static void eigenfunction(arb_t res, arb_t z, arb_t phi,
                                       arb_ptr coefs, slong N, arb_t nu,
                                       arb_t mu0, int (*index)(int),
                                       slong prec) {
  arb_t mu, term, tmp;
  slong prec_local;

  arb_init(mu);
  arb_init(term);
  arb_init(tmp);

  for (slong i = 0; i < N; i++) {
    arb_mul_si(mu, mu0, index(i), prec);

    prec_local = prec;
    do {
      arb_hypgeom_legendre_p(term, nu, mu, z, 0, prec_local);
      prec_local *=2;
    } while (!arb_is_finite(term) && prec_local <= 16*prec);

    arb_mul(tmp, mu, phi, prec);
    arb_sin(tmp, tmp, prec);

    arb_mul(term, term, tmp, prec);
    arb_mul(term, term, coefs + i, prec);

    arb_add(res, res, term, prec);
  }

  arb_clear(mu);
  arb_clear(term);
  arb_clear(tmp);
}

static void maximize(arb_t max, arb_ptr coefs, slong N, arb_ptr v1,
                     arb_ptr v2, arb_t critical_point, arb_t nu, arb_t mu0,
                     int (*index)(int), slong prec) {
  arb_ptr evals;
  arb_t z, phi;
  arf_t *intervals_low, *next_intervals_low, *intervals_upp, *next_intervals_upp;
  arf_t t_low, t_upp, max_low, max_upp, tmp;
  mag_t tol;
  slong num_intervals, next_num_intervals;
  int k, done;

  intervals_low = new arf_t[1];
  intervals_upp = new arf_t[1];
  arf_init(intervals_low[0]);
  arf_init(intervals_upp[0]);

  arb_init(z);
  arb_init(phi);

  arf_init(t_low);
  arf_init(t_upp);
  arf_init(max_low);
  arf_init(max_upp);
  arf_init(tmp);

  mag_init(tol);

  mag_set_ui_2exp_si(tol, 1, -16);

  arf_set_si(intervals_low[0], 0);
  arf_set_si(intervals_upp[0], 1);

  num_intervals = 1;

  k = 0;
  done = 0;

  while (!done) {
    evals = _arb_vec_init(num_intervals);
    next_intervals_low = new arf_t[num_intervals*2];
    next_intervals_upp = new arf_t[num_intervals*2];

    arf_zero(max_low);
    arf_zero(max_upp);

    k++;
    next_num_intervals = 0;

    /* Compute the eigenfunction on the intervals and update the bound
     * for the maximum */
    for (slong i = 0; i < num_intervals; i++) {
      parametrization(z, phi, intervals_low[i], intervals_upp[i], v1, v2,
                       critical_point, prec);
      eigenfunction(evals + i, z, phi, coefs, N, nu, mu0, index, prec);
      if (arb_is_finite(evals + i)){
        arb_get_abs_ubound_arf(tmp, evals + i, prec);
        arf_max(max_upp, max_upp, tmp);
        arb_get_abs_lbound_arf(tmp, evals + i, prec);
        arf_max(max_low, max_low, tmp);
      }
    }

    /* Check if done */
    arf_sub(tmp, max_upp, max_low, prec, ARF_RND_UP);
    if (!(k < 15) ||
        ((arf_cmpabs_mag(tmp, tol) <= 0) && arf_cmp_si(max_upp, 0) > 0)) {
      done = 1;
    } else {
      /* Split the intervals where the maximum could be obtained */
      for (slong i = 0; i < num_intervals; i++) {
        arb_get_abs_ubound_arf(tmp, evals + i, prec);
        if (arf_is_nan(tmp) || (arf_cmp(tmp, max_low) >= 0)) {
          /* Split the interval */
          arf_init(next_intervals_low[next_num_intervals]);
          arf_init(next_intervals_low[next_num_intervals + 1]);
          arf_init(next_intervals_upp[next_num_intervals]);
          arf_init(next_intervals_upp[next_num_intervals + 1]);

          arf_add(tmp, intervals_low[i], intervals_upp[i], prec, ARF_RND_NEAR);
          arf_div_si(tmp, tmp, 2, prec, ARF_RND_NEAR);
          arf_set(next_intervals_low[next_num_intervals], intervals_low[i]); // Lower
          arf_set(next_intervals_upp[next_num_intervals], tmp);              // Mid
          arf_set(next_intervals_low[next_num_intervals + 1], tmp);          // Mid
          arf_set(next_intervals_upp[next_num_intervals + 1], intervals_upp[i]); // Upper

          next_num_intervals += 2;
        }
      }
    }

    _arb_vec_clear(evals, num_intervals);

    for (slong i = 0; i < num_intervals; i++) {
      arf_clear(intervals_low[i]);
      arf_clear(intervals_upp[i]);
    }
    delete [] intervals_low;
    delete [] intervals_upp;

    intervals_low = next_intervals_low;
    intervals_upp = next_intervals_upp;

    num_intervals = next_num_intervals;
  }

  arb_set_arf(max, max_upp);

  for (slong i = 0; i < num_intervals; i++) {
    arf_clear(intervals_low[i]);
    arf_clear(intervals_upp[i]);
  }

  delete [] intervals_low;
  delete [] intervals_upp;

  arb_clear(z);
  arb_clear(phi);

  arf_clear(t_low);
  arf_clear(t_upp);
  arf_clear(max_low);
  arf_clear(max_upp);
  arf_clear(tmp);
}

void enclose(mpfr_t eps_mpfr, int angles_coefs[], mpfr_t *coefs_mpfr, int N,
             mpfr_t nu_mpfr, int (*index)(int)) {
  arb_ptr v1, v2, coefs;
  arb_t eps, nu, mu0, theta_bound, critical_point, norm, max, tmp;
  arf_t eps_upp;
  slong prec;

  v1 = _arb_vec_init(3);
  v2 = _arb_vec_init(3);
  coefs = _arb_vec_init(N);

  arb_init(eps);
  arb_init(nu);
  arb_init(mu0);
  arb_init(theta_bound);
  arb_init(critical_point);
  arb_init(norm);
  arb_init(max);
  arb_init(tmp);

  arf_init(eps_upp);

  prec = mpfr_get_default_prec();

  for (slong i = 0; i < N; i++) {
    arf_set_mpfr(arb_midref(coefs + i), coefs_mpfr[i]);
  }

  arf_set_mpfr(arb_midref(nu), nu_mpfr);

  /* Compute enclosures of the required parameters */
  angles_to_vectors_arb(v1, v2, theta_bound, critical_point, angles_coefs, prec);

  arb_set_si(mu0, -angles_coefs[1]);
  arb_div_si(mu0, mu0, angles_coefs[0], prec);

  /* Compute an enclosure of a lower bound of the norm */
  integral_norm(norm, coefs, N, theta_bound, nu, mu0, index, prec);

  maximize(max, coefs, N, v1, v2, critical_point, nu, mu0, index, prec);

  /* Set eps equal to the square root of the area of the triangle */
  for (slong i = 0; i < 3; i++) {
    arb_set_si(tmp, angles_coefs[2*i]);
    arb_div_si(tmp, tmp, angles_coefs[2*i + 1], prec);
    arb_add(eps, eps, tmp, prec);
  }
  arb_const_pi(tmp, prec);
  arb_sub_si(eps, eps, 1, prec);
  arb_mul(eps, eps, tmp, prec);
  arb_sqrt(eps, eps, prec);

  /* Multiply with the maximum and divide with the norm */
  arb_mul(eps, eps, max, prec);
  arb_div(eps, eps, norm, prec);

  arb_get_abs_ubound_arf(eps_upp, eps, prec);
  arf_get_mpfr(eps_mpfr, eps_upp, MPFR_RNDU);

  _arb_vec_clear(v1, 3);
  _arb_vec_clear(v2, 3);
  _arb_vec_clear(coefs, N);
  arb_clear(eps);
  arb_clear(nu);
  arb_clear(mu0);
  arb_clear(theta_bound);
  arb_clear(critical_point);
  arb_clear(norm);
  arb_clear(max);
  arb_clear(tmp);

  arf_clear(eps_upp);
}
