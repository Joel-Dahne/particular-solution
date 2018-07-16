#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_calc.h"

static void
angles_to_vectors_arb(arb_ptr v1, arb_ptr v2, arb_t theta_bound_low,
                                  arb_t theta_bound_upp, arb_t critical_point,
                                  int angles_coefs[], slong prec)
{
  arb_ptr angles, w;
  arb_t S, theta1, theta2, tmp, tmp_dot;

  angles = _arb_vec_init(3);
  w = _arb_vec_init(3);
  arb_init(S);
  arb_init(theta1);
  arb_init(theta2);
  arb_init(tmp);
  arb_init(tmp_dot);

  for (slong i = 0; i < 3; i++)
  {
    arb_set_si(angles + i, angles_coefs[2*i]);
    arb_div_si(angles + i, angles + i, angles_coefs[2*i + 1], prec);
  }

  /* Compute the sum of the angles */
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

  /* Compute theta2 for v2 */
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
  /* w_3*dot(v, w) - v_3*dot(w, w)*/
  arb_submul(tmp, v1 + 2, tmp_dot, prec);
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


  /* Compute the minimum theta value */
  arb_min(theta_bound_low, theta1, theta2, prec);
  arb_min(theta_bound_low, theta_bound_low, tmp, prec);
  /* Compute the maximum theta value */
  arb_max(theta_bound_upp, theta1, theta2, prec);
  arb_max(theta_bound_upp, theta_bound_upp, tmp, prec);

  _arb_vec_clear(angles, 3);
  _arb_vec_clear(w, 3);
  arb_clear(S);
  arb_clear(theta1);
  arb_clear(theta2);
  arb_clear(tmp);
  arb_clear(tmp_dot);
}

static int
integrand(acb_ptr out, const acb_t inp, void *param_void, slong order,
                     slong prec)
{
  acb_ptr param;
  acb_t tmp;
  slong prec_local;

  acb_init(tmp);

  /* param = [nu, mu, coef] */
  param = (acb_ptr)(param_void);

  if (order > 1)
    flint_abort();

  acb_cos(out, inp, prec);

  if (order == 1)
  {
    acb_sub_si(tmp, out, 1, prec);
    if (!arb_is_nonpositive(acb_realref(tmp)))
    {
      acb_indeterminate(out);
      return 0;
    }
  }

  prec_local = prec;
  do
  {
    acb_hypgeom_legendre_p(tmp, param, param + 1, out, 0, prec_local);
    prec_local *= 2;
  } while (!acb_is_finite(tmp) && prec_local <= 16*prec);

  acb_sqr(out, tmp, prec);
  acb_sin(tmp, inp, prec);
  acb_mul(out, out, tmp, prec);
  acb_mul(out, out, param + 2, prec);

  acb_clear(tmp);

  return 0;
}

static void
integral_norm(arb_t norm, arb_ptr coefs, int N, arb_t theta_bound,
                          arb_t nu, arb_t mu0, int (*index)(int), slong prec)
{
  acb_t zero, theta_bound_c, norm_part;
  acb_ptr params;
  arb_t integral_phi, mu;
  mag_t tol;
  slong prec_local;

  params = _acb_vec_init(3);

  acb_init(zero);
  acb_init(theta_bound_c);
  acb_init(norm_part);
  arb_init(integral_phi);
  arb_init(mu);

  mag_init(tol);

  arb_set(acb_realref(params + 0), nu); /* params[0] = nu */

  arb_set_str(acb_realref(zero), "1e-1", prec);
  arb_set(acb_realref(theta_bound_c), theta_bound);

  /* Compute the integral corresponding to phi */
  arb_const_pi(integral_phi, prec);
  arb_div_si(integral_phi, integral_phi, 2, prec);
  arb_div(integral_phi, integral_phi, mu0, prec);
  arb_neg(integral_phi, integral_phi);

  for (slong i = 0; i < N; i++)
  {
    /* params[1] = i*mu0 */
    arb_mul_si(acb_realref(params + 1), mu0, index(i), prec);
    /* params[2] = sqr(coefs[i]) */
    arb_sqr(acb_realref(params + 2), coefs + i, prec);

    /* Compute the integral corresponding to theta */
    prec_local = prec;
    do
    {
      acb_calc_integrate(norm_part, integrand, params, zero, theta_bound_c,
                         prec, tol, NULL, prec_local);
      prec_local *= 2;
    } while (!acb_is_finite(norm_part));

    /* Add to the norm */
    arb_add(norm, norm, acb_realref(norm_part), prec);
  }
  /* Multiply with the integral corresponding to phi */
  arb_mul(norm, norm, integral_phi, prec);

  arb_sqrt(norm, norm, prec);

  _acb_vec_clear(params, 3);

  acb_clear(zero);
  acb_clear(theta_bound_c);
  acb_clear(norm_part);
  arb_clear(integral_phi);
  arb_clear(mu);
}

static void
parametrization(arb_ptr z, arb_ptr phi, arb_t t, arb_ptr v1, arb_ptr v2,
                slong order, slong prec)
{
  arb_ptr w, x, y, dx, dy, tmp1, tmp2;

  w = _arb_vec_init(3);
  x = _arb_vec_init(order);
  y = _arb_vec_init(order);
  dx = _arb_vec_init(order - 1);
  dy = _arb_vec_init(order - 1);
  tmp1 = _arb_vec_init(order);
  tmp2 = _arb_vec_init(order);

  _arb_vec_zero(z, order);

  /* w = v2 - v1 */
  _arb_vec_sub(w, v2, v1, 3, prec);

  /* x = v1_x + t*w_x */
  /* x' = w_x */
  arb_mul(x, t, w, prec);
  arb_add(x, x, v1, prec);
  if (order > 1)
    arb_set(x + 1, w);

  /* y = v1_y + t*w_y */
  /* y' = w_y */
  arb_mul(y, t, w + 1, prec);
  arb_add(y, y, v1 + 1, prec);
  if (order > 1)
    arb_set(y + 1, w + 1);

  /* z = v1_z + t*w_z */
  /* z' = w_z */
  arb_mul(z, t, w + 2, prec);
  arb_add(z, z, v1 + 2, prec);
  if (order > 1)
    arb_set(z + 1, w + 2);

  /* tmp1 = x^2 */
  _arb_poly_mullow(tmp1, x, order, x, order, order, prec);
  /* tmp2 = y^2 */
  _arb_poly_mullow(tmp2, y, order, y, order, order, prec);
  /* tmp1 = x^2 + y^2 */
  _arb_vec_add(tmp1, tmp1, tmp2, order, prec);
  /* tmp2 = z^2 */
  _arb_poly_mullow(tmp2, z, order, z, order, order, prec);
  /* tmp1 = x^2 + y^2 + z^2 */
  _arb_vec_add(tmp1, tmp1, tmp2, order, prec);
  /* tmp2 = sqrt(x^2 + y^2 + z^2) */
  _arb_poly_sqrt_series(tmp2, tmp1, order, order, prec);

  /* x, y, z = z, y, z / ||(x, y, z)|| */
  _arb_poly_div_series(tmp1, x, order, tmp2, order, order, prec);
  _arb_vec_swap(tmp1, x, order);
  _arb_poly_div_series(tmp1, y, order, tmp2, order, order, prec);
  _arb_vec_swap(tmp1, y, order);
  _arb_poly_div_series(tmp1, z, order, tmp2, order, order, prec);
  _arb_vec_swap(tmp1, z, order);

  /* Compute the Taylor expansion of phi' */
  _arb_poly_derivative(dx, x, order, prec);
  _arb_poly_derivative(dy, y, order, prec);

  _arb_poly_mullow(tmp1, x, order - 1, dy, order - 1, order - 1, prec);
  _arb_poly_mullow(tmp2, y, order - 1, dx, order - 1, order - 1, prec);
  _arb_vec_sub(phi, tmp1, tmp2, order - 1, prec);
  _arb_poly_mullow(tmp1, x, order - 1, x, order - 1, order - 1, prec);
  _arb_poly_mullow(tmp2, y, order - 1, y, order - 1, order - 1, prec);
  _arb_vec_add(tmp1, tmp1, tmp2, order - 1, prec);
  _arb_poly_div_series(tmp2, phi, order - 1, tmp1, order - 1, order - 1, prec);
  _arb_vec_swap(phi, tmp2, order - 1);
  _arb_poly_integral(phi, phi, order, prec);

  /* phi = atan2(y, x) */
  arb_atan2(phi, y, x, prec);

  _arb_vec_clear(w, 3);
  _arb_vec_clear(x, order);
  _arb_vec_clear(y, order);
  _arb_vec_clear(dx, order - 1);
  _arb_vec_clear(dy, order - 1);
  _arb_vec_clear(tmp1, order);
  _arb_vec_clear(tmp2, order);
}

static void
legendre_taylor(arb_ptr res, arb_ptr z, slong order, arb_t nu, arb_t mu,
                slong prec)
{
  arb_ptr legendre, tmp_poly;
  arb_t z2, nu2, mu2, onemz2, tmp1, tmp2;
  slong prec_local, n;

  legendre = _arb_vec_init(order);
  tmp_poly = _arb_vec_init(order);

  arb_init(z2);
  arb_init(nu2);
  arb_init(mu2);
  arb_init(onemz2);
  arb_init(tmp1);
  arb_init(tmp2);

  arb_sqr(z2, z, prec);
  arb_sqr(nu2, nu, prec);
  arb_sqr(mu2, mu, prec);
  arb_neg(onemz2, z2);
  arb_add_si(onemz2, onemz2, 1, prec);

  if (order > 0)
  {
    /* a_0 */
    prec_local = prec;
    do {
      arb_hypgeom_legendre_p(legendre, nu, mu, z, 0, prec_local);
      prec_local = 2*prec_local;
    } while (!arb_is_finite(legendre) && prec_local < 16*prec);

    if (order > 1)
    {
      /* a_1 */
      arb_add_si(tmp1, nu, 1, prec);
      do {
        arb_hypgeom_legendre_p(legendre + 1, tmp1, mu, z, 0, prec_local);
        prec_local = 2*prec_local;
      } while (!arb_is_finite(legendre + 1) && prec_local < 16*prec);
      arb_sub(tmp1, tmp1, mu, prec);
      arb_mul(legendre + 1, legendre + 1, tmp1, prec);

      arb_add_si(tmp1, nu, 1, prec);
      arb_mul(tmp1, tmp1, z, prec);
      arb_mul(tmp1, tmp1, legendre, prec);
      arb_sub(legendre + 1, tmp1, legendre + 1, prec);

      arb_div(legendre + 1, legendre + 1, onemz2, prec);

      if (order > 2)
      {
        /* a_2 */
        arb_mul_si(tmp1, z, 2, prec);
        arb_mul(legendre + 2, tmp1, legendre + 1, prec);

        arb_div(tmp1, mu2, onemz2, prec);
        arb_add_si(tmp2, nu, 1, prec);
        arb_mul(tmp2, tmp2, nu, prec);
        arb_sub(tmp1, tmp2, tmp1, prec);
        arb_submul(legendre + 2, tmp1, legendre, prec);

        arb_div(legendre + 2, legendre + 2, onemz2, prec);
        arb_div_si(legendre + 2, legendre + 2, 2, prec);

        if (order > 3)
        {
          /* a_3 */
          arb_add(tmp1, nu2, nu, prec);
          arb_mul(tmp1, tmp1, z2, prec);
          arb_mul_si(tmp2, mu2, 2, prec);
          arb_add(tmp1, tmp1, tmp2, prec);
          arb_sub(tmp1, tmp1, nu2, prec);
          arb_sub(tmp1, tmp1, nu, prec);
          arb_mul(tmp1, tmp1, z, prec);
          arb_mul_si(tmp1, tmp1, 2, prec);
          arb_mul(legendre + 3, tmp1, legendre + 0, prec);

          arb_add(tmp1, nu2, nu, prec);
          arb_add_si(tmp1, tmp1, 2, prec);
          arb_mul(tmp1, tmp1, z2, prec);
          arb_add(tmp1, tmp1, mu2, prec);
          arb_sub(tmp1, tmp1, nu2, prec);
          arb_sub(tmp1, tmp1, nu, prec);
          arb_add_si(tmp1, tmp1, 2, prec);
          arb_add_si(tmp2, z, 1, prec);
          arb_mul(tmp1, tmp1, tmp2, prec);
          arb_add_si(tmp2, z, -1, prec);
          arb_mul(tmp1, tmp1, tmp2, prec);
          arb_submul(legendre + 3, tmp1, legendre + 1, prec);

          arb_sqr(tmp1, z2, prec);
          arb_mul_si(tmp2, z2, -4, prec);
          arb_mul_si(tmp1, tmp1, 2, prec);
          arb_add(tmp1, tmp1, tmp2, prec);
          arb_add_si(tmp1, tmp1, 2, prec);
          arb_mul(tmp1, tmp1, z, prec);
          arb_mul_si(tmp1, tmp1, 2, prec);
          arb_addmul(legendre + 3, tmp1, legendre + 2, prec);

          arb_pow_ui(tmp1, onemz2, 3, prec);
          arb_mul_si(tmp1, tmp1, 6, prec);
          arb_div(legendre + 3, legendre + 3, tmp1, prec);

          /* Compute the remaining terms using the recurrence relation */
          for (slong i = 4; i < order; i++)
          {
            n = i - 4;

            /* Term for n */
            arb_add_si(tmp1, nu, n + 1, prec);
            arb_add_si(tmp2, nu, -n, prec);
            arb_mul(tmp1, tmp1, tmp2, prec);
            arb_mul(legendre + n + 4, tmp1, legendre + n, prec);

            /* Term for n + 1 */
            arb_add(tmp1, nu2, nu, prec);
            arb_neg(tmp1, tmp1);
            arb_add_si(tmp1, tmp1, 5*n + 3, prec);
            arb_div_si(tmp1, tmp1, 2, prec);
            arb_add_si(tmp1, tmp1, n*n, prec);
            arb_mul_si(tmp1, tmp1, -4, prec);
            arb_mul(tmp1, tmp1, z, prec);
            arb_addmul(legendre + n + 4, tmp1, legendre + n + 1, prec);

            /* Term for n + 2 */
            arb_add(tmp1, nu, nu2, prec);
            arb_add_si(tmp1, tmp1, -6*n*n - 24*n - 24, prec);
            arb_mul(tmp1, tmp1, z2, prec);
            arb_add(tmp1, tmp1, mu2, prec);
            arb_sub(tmp1, tmp1, nu2, prec);
            arb_sub(tmp1, tmp1, nu, prec);
            arb_add_si(tmp1, tmp1, 2*n*n + 8*n + 8, prec);
            arb_addmul(legendre + n + 4, tmp1, legendre + n + 2, prec);

            /* Term for n + 3 */
            arb_add_si(tmp1, z, 1, prec);
            arb_mul(tmp1, tmp1, z, prec);
            arb_add_si(tmp2, z, - 1, prec);
            arb_mul(tmp1, tmp1, tmp2, prec);
            arb_mul_si(tmp1, tmp1, -2*(n + 3)*(2*n + 5), prec);
            arb_addmul(legendre + n + 4, tmp1, legendre + n + 3, prec);

            /* Factor */
            arb_add_si(tmp1, z, -1, prec);
            arb_sqr(tmp1, tmp1, prec);
            arb_add_si(tmp2, z, 1, prec);
            arb_sqr(tmp2, tmp2, prec);
            arb_mul(tmp1, tmp1, tmp2, prec);
            arb_mul_si(tmp1, tmp1, (n + 4)*(n + 3), prec);
            arb_div(legendre + n + 4, legendre + n + 4, tmp1, prec);
          }
        }
      }
    }
  }

  /* Compose the Taylor series for the Legendre function with that of
   * z */
  _arb_vec_set(tmp_poly, z, order);
  arb_zero(tmp_poly);
  _arb_poly_compose_series(res, legendre, order, tmp_poly, order, order, prec);

  _arb_vec_clear(legendre, order);
  _arb_vec_clear(tmp_poly, order);

  arb_clear(z2);
  arb_clear(nu2);
  arb_clear(mu2);
  arb_clear(onemz2);
  arb_clear(tmp1);
  arb_clear(tmp2);
}

static void
eigenfunction_taylor(arb_ptr res, arb_ptr z, arb_ptr phi, arb_ptr coefs,
                     slong order, slong N, arb_t nu, arb_t mu0, int (*index)(int),
                     slong prec)
{
  arb_ptr term_sin, term_legendre, tmp;
  arb_t mu;

  term_sin = _arb_vec_init(order);
  term_legendre = _arb_vec_init(order);
  tmp = _arb_vec_init(order);

  arb_init(mu);

  _arb_vec_zero(res, order);

  for (slong i = 0; i < N; i++)
  {
    arb_mul_si(mu, mu0, index(i), prec);
    _arb_vec_scalar_mul(tmp, phi, order, mu, prec);
    _arb_poly_sin_series(term_sin, tmp, order, order, prec);
    legendre_taylor(term_legendre, z, order, nu, mu, prec);
    _arb_poly_mullow(tmp, term_sin, order, term_legendre, order, order, prec);
    _arb_vec_scalar_addmul(res, tmp, order, coefs + i, prec);
  }

  _arb_vec_clear(term_sin, order);
  _arb_vec_clear(term_legendre, order);
  _arb_vec_clear(tmp, order);

  arb_clear(mu);
}

static void
maximize_taylor(arb_t max, arb_t t, arb_ptr coefs, slong N, slong order,
                arb_ptr v1, arb_ptr v2, arb_t nu, arb_t mu0, int (*index)(int),
                slong prec) {
  arb_ptr z, phi, poly;
  arb_t x, rest_term;

  z = _arb_vec_init(order + 1);
  phi = _arb_vec_init(order + 1);
  poly = _arb_vec_init(order + 1);

  arb_init(x);
  arb_init(rest_term);

  /* Compute the Taylor polynomial of the eigenfunction at the midpoint of t */
  arb_set_arf(x, arb_midref(t));
  parametrization(z, phi, x, v1, v2, order, prec);
  eigenfunction_taylor(poly, z, phi, coefs, order, N, nu, mu0,
                       index, prec);

  /* Enclose the Taylor polynomial evaluated at t - mid(t) */
  arb_sub(x, t, x, prec);
  _arb_poly_evaluate(max, poly, order, x, prec);
  arb_abs(max, max);

  /* Compute the rest term of the Taylor expansion */
  parametrization(z, phi, t, v1, v2, order + 1, prec);
  eigenfunction_taylor(poly, z, phi, coefs, order + 1, N, nu, mu0,
                       index, prec);

  arb_pow_ui(rest_term, x, order, prec);
  arb_mul(rest_term, rest_term, poly + order, prec);

  /* Add the rest term to the maximum */
  arb_add(max, max, rest_term, prec);

  _arb_vec_clear(z, order + 1);
  _arb_vec_clear(phi, order + 1);
  _arb_vec_clear(poly, order);

  arb_clear(x);
  arb_clear(rest_term);
}

static void
maximize(arb_t max, arb_ptr coefs, slong N, arb_ptr v1, arb_ptr v2,
         arb_t nu, arb_t mu0, int (*index)(int), slong prec) {
  arb_ptr evals;
  arf_t *intervals_low, *next_intervals_low, *intervals_upp, *next_intervals_upp;
  arb_t t;
  arf_t t_low, t_upp, max_low, max_upp, tmp;
  slong order, num_intervals, next_num_intervals;
  int k, done;

  order = 16;

  intervals_low = new arf_t[1];
  intervals_upp = new arf_t[1];
  arf_init(intervals_low[0]);
  arf_init(intervals_upp[0]);

  arb_init(t);

  arf_init(t_low);
  arf_init(t_upp);
  arf_init(max_low);
  arf_init(max_upp);
  arf_init(tmp);

  arf_set_d(intervals_low[0], 0);
  arf_set_d(intervals_upp[0], 1);

  num_intervals = 1;

  k = 0;
  done = 0;

  while (!done)
  {
    evals = _arb_vec_init(num_intervals);
    next_intervals_low = new arf_t[num_intervals*2];
    next_intervals_upp = new arf_t[num_intervals*2];

    k++;
    next_num_intervals = 0;

    arb_zero(max);
    arf_zero(max_upp);
    arf_zero(max_low);

    /* Compute a bound for the maximum for the eigenfunction on the
     * intervals */
    for (slong i = 0; i < num_intervals; i++)
    {
      arb_set_interval_arf(t, intervals_low[i], intervals_upp[i], prec);

      maximize_taylor(evals + i, t, coefs, N, order, v1, v2, nu, mu0, index,
                      prec);

      arb_get_abs_ubound_arf(tmp, evals + i, prec);
      arf_max(max_upp, max_upp, tmp);

      if (arb_is_finite(evals + i))
      {
        arb_get_abs_lbound_arf(tmp, evals + i, prec);
        arf_max(max_low, max_low, tmp);
      }
    }

    arb_set_interval_arf(max, max_low, max_upp, prec);

    /* Check if done */
    if (!(k < 30) || (arb_rel_error_bits(max) < -5))
    {
      done = 1;
    }
    else
    {
      /* Split the intervals where the maximum could be obtained */
      for (slong i = 0; i < num_intervals; i++)
      {
        arb_get_abs_ubound_arf(tmp, evals + i, prec);
        if (!(arf_cmp(tmp, max_low) < 0))
        {
          /* Split the interval */
          arf_init(next_intervals_low[next_num_intervals]);
          arf_init(next_intervals_low[next_num_intervals + 1]);
          arf_init(next_intervals_upp[next_num_intervals]);
          arf_init(next_intervals_upp[next_num_intervals + 1]);

          arf_add(tmp, intervals_low[i], intervals_upp[i], prec, ARF_RND_NEAR);
          arf_div_si(tmp, tmp, 2, prec, ARF_RND_NEAR);

          arf_set(next_intervals_low[next_num_intervals], intervals_low[i]);
          arf_set(next_intervals_upp[next_num_intervals], tmp);
          arf_set(next_intervals_low[next_num_intervals + 1], tmp);
          arf_set(next_intervals_upp[next_num_intervals + 1], intervals_upp[i]);

          next_num_intervals += 2;
        }
      }
    }

    _arb_vec_clear(evals, num_intervals);

    for (slong i = 0; i < num_intervals; i++)
    {
      arf_clear(intervals_low[i]);
      arf_clear(intervals_upp[i]);
    }
    delete [] intervals_low;
    delete [] intervals_upp;

    intervals_low = next_intervals_low;
    intervals_upp = next_intervals_upp;

    num_intervals = next_num_intervals;
  }

  for (slong i = 0; i < num_intervals; i++)
  {
    arf_clear(intervals_low[i]);
    arf_clear(intervals_upp[i]);
  }

  delete [] intervals_low;
  delete [] intervals_upp;

  arb_clear(t);

  arf_clear(t_low);
  arf_clear(t_upp);
  arf_clear(max_low);
  arf_clear(max_upp);
  arf_clear(tmp);
}

void
enclose(mpfr_t eps_mpfr, int angles_coefs[], mpfr_t *coefs_mpfr, int N,
        mpfr_t nu_mpfr, int (*index)(int)) {
  arb_ptr v1, v2, coefs;
  arb_t eps, nu, mu0, theta_bound_low, theta_bound_upp, critical_point, norm,
    max, eigenvalue, tmp;
  arf_t eps_upp;
  slong prec;

  v1 = _arb_vec_init(3);
  v2 = _arb_vec_init(3);
  coefs = _arb_vec_init(N);

  arb_init(eps);
  arb_init(nu);
  arb_init(mu0);
  arb_init(theta_bound_low);
  arb_init(theta_bound_upp);
  arb_init(critical_point);
  arb_init(norm);
  arb_init(max);
  arb_init(eigenvalue);
  arb_init(tmp);

  arf_init(eps_upp);

  prec = mpfr_get_default_prec();

  for (slong i = 0; i < N; i++)
  {
    arf_set_mpfr(arb_midref(coefs + i), coefs_mpfr[i]);
  }

  arf_set_mpfr(arb_midref(nu), nu_mpfr);

  /* Compute enclosures of the required parameters */
  angles_to_vectors_arb(v1, v2, theta_bound_low, theta_bound_upp,
                        critical_point, angles_coefs, prec);

  arb_set_si(mu0, -angles_coefs[1]);
  arb_div_si(mu0, mu0, angles_coefs[0], prec);

  /* Compute an enclosure of a lower bound of the norm */
  integral_norm(norm, coefs, N, theta_bound_low, nu, mu0, index, prec);
  maximize(max, coefs, N, v1, v2, nu, mu0, index, prec);

  /* Set eps equal to the square root of the area of the triangle */
  for (slong i = 0; i < 3; i++)
  {
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

  /* Compute lower and upper bounds for the eigenvalue */
  arb_add_si(eigenvalue, nu, 1, prec);
  arb_mul(eigenvalue, eigenvalue, nu, prec);
  arb_set_si(tmp, 1);
  arb_add_error(tmp, eps);
  arb_div(eigenvalue, eigenvalue, tmp, prec);

  flint_printf(" ");
  arb_printn(eigenvalue, (slong)ceil(prec*log10(2)), 0);

  arb_get_abs_ubound_arf(eps_upp, eps, prec);
  arf_get_mpfr(eps_mpfr, eps_upp, MPFR_RNDU);

  _arb_vec_clear(v1, 3);
  _arb_vec_clear(v2, 3);
  _arb_vec_clear(coefs, N);
  arb_clear(eps);
  arb_clear(nu);
  arb_clear(mu0);
  arb_clear(theta_bound_low);
  arb_clear(theta_bound_upp);
  arb_clear(critical_point);
  arb_clear(norm);
  arb_clear(max);
  arb_clear(eigenvalue);
  arb_clear(tmp);

  arf_clear(eps_upp);
}
