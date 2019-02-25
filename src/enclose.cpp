#include "enclose.h"

#include "norm.h"
#include "eigenfunction.h"
#include "arb_poly.h"

static void
maximize_taylor(arb_t max, geom_t geom, arb_t t, arb_ptr coefs, slong N,
                slong order, arb_t nu, slong prec) {
  arb_ptr z, phi, poly;
  arb_t x, rest_term;

  z = _arb_vec_init(order + 1);
  phi = _arb_vec_init(order + 1);
  poly = _arb_vec_init(order + 1);

  arb_init(x);
  arb_init(rest_term);

  /* Compute the Taylor polynomial of the eigenfunction at the midpoint of t */
  arb_set_arf(x, arb_midref(t));
  parametrization(z, phi, geom, x, order, 0, prec);
  eigenfunction_series(poly, geom, coefs, N, z, phi, nu, 0, order, prec);

  /* Enclose the Taylor polynomial evaluated at t - mid(t) */
  arb_sub(x, t, x, prec);
  _arb_poly_evaluate(max, poly, order, x, prec);
  arb_abs(max, max);

  /* Compute the rest term of the Taylor expansion */
  parametrization(z, phi, geom, t, order + 1, 0, prec);
  eigenfunction_series(poly, geom, coefs, N, z, phi, nu, 0, order + 1, prec);

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
maximize(arb_t max, geom_t geom, arb_ptr coefs, slong N, arb_t nu, slong prec) {
  arb_ptr evals;
  arf_t *intervals_low, *next_intervals_low, *intervals_upp, *next_intervals_upp;
  arb_t t;
  arf_t t_low, t_upp, max_low, max_upp, tmp;
  slong order, num_intervals, next_num_intervals, total_intervals, limit;
  int k, done;

  if (N > 2)
  {
    order = N;
  }
  else
  {
    order = 2;
  }

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
  total_intervals = num_intervals;
  limit = 500;

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

      maximize_taylor(evals + i, geom, t, coefs, N, order, nu, prec);

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
    total_intervals += num_intervals;

    if (total_intervals > limit)
    {
      prec *= 2;
      limit *= 2;
    }
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
enclose(arb_t nu_enclosure, geom_t geom, arb_ptr coefs,
        slong N, arb_t nu, slong prec) {
  arb_t eps, norm, max, eigenvalue, tmp;

  arb_init(eps);
  arb_init(norm);
  arb_init(max);
  arb_init(eigenvalue);
  arb_init(tmp);

  /* Compute an enclosure of a lower bound of the norm */
  integral_norm(norm, geom, coefs, N, nu, 0, prec);

  /* Compute an enclosure of the maximum on the boundary */
  maximize(max, geom, coefs, N, nu, prec);

  /* Set eps equal to the square root of the area of the triangle */
  for (slong i = 0; i < 3; i++)
  {
    arb_set_fmpq(tmp, geom->angles + i, prec);
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

  /* Compute lower and upper bounds for the enclosure of nu */
  arb_set_d(tmp, 0.25);
  arb_add(tmp, tmp, eigenvalue, prec);
  arb_sqrt(tmp, tmp, prec);
  arb_set_d(eps, 0.5);
  arb_sub(tmp, tmp, eps, prec);

  if (arb_is_finite(tmp) && arb_overlaps(nu_enclosure, tmp))
  {
    arb_intersection(nu_enclosure, nu_enclosure, tmp, prec);
  }

  arb_clear(eps);
  arb_clear(norm);
  arb_clear(max);
  arb_clear(eigenvalue);
  arb_clear(tmp);
}
