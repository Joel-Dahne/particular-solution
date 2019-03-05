#include "maximum.h"

#include "eigenfunction.h"
#include "arb_poly.h"

/* Enclose the maximum value of the eigenfunction with by the given
 * coefficients and parameters on the interval t. The argument vertex
 * indicates which vertex the eigenfunction should be computed from, t
 * is then used when computing the parameterization of the opposite
 * boundary.
 *
 * The computations are done by computing a Taylor expansion
 * with n terms centered at the midpoint of the interval t for which
 * the maximum, including the enclosed error term, is computed.
 *
 * TODO: At the moment we compute the maximum of the enclosure of the
 * Taylor expansion on the while interval and the enclosure around the
 * midpoint. We do not try to find the maximum value. While this is
 * correct, i.e. the enclosure it gives contains the maximum on the
 * interval, it is not optimal. It would be better if we were able to
 * enclose only the maximum and do so in an efficient way. */
void
maximize_series(arb_t max, geom_t geom, arb_t t, arb_ptr coefs, slong N,
                arb_t nu, slong n, slong vertex, slong prec)
{
  arb_ptr z, phi, series;
  arb_t enclosure, rest_term, tmp;

  z = _arb_vec_init(n + 1);
  phi = _arb_vec_init(n + 1);
  series = _arb_vec_init(n + 1);

  arb_init(enclosure);
  arb_init(rest_term);
  arb_init(tmp);

  /* Compute the Taylor polynomial of the eigenfunction at the midpoint of t */
  /* tmp = mid(t) */
  arb_set_arf(tmp, arb_midref(t));
  parametrization(z, phi, geom, tmp, n, vertex, prec);
  eigenfunction_series(series, geom, coefs, N, z, phi, nu, vertex, n, prec);

  /* The max is lower bounded by the value at the midpoint. */
  arb_abs(max, series + 0);

  /* Enclose the Taylor polynomial evaluated at t - mid(t) */
  /* tmp = t - mid(t) */
  arb_sub(tmp, t, tmp, prec);
  _arb_poly_evaluate(tmp, series, n, tmp, prec);
  arb_abs(enclosure, tmp);

  /* Compute the rest term of the Taylor expansion */
  parametrization(z, phi, geom, t, n + 1, vertex, prec);
  eigenfunction_series(series, geom, coefs, N, z, phi, nu, vertex, n + 1, prec);

  /* rest_term = (t - mid(t))^n*(enclosure of coefficient n) */
  arb_pow_ui(rest_term, tmp, n, prec);
  arb_mul(rest_term, rest_term, series + n, prec);

    /* Add the rest term to the maximum */
  arb_add(enclosure, enclosure, rest_term, prec);

  /* Get the maximum of the midpoint and the enclosure. */
  arb_max(max, max, enclosure, prec);

  _arb_vec_clear(z, n + 1);
  _arb_vec_clear(phi, n + 1);
  _arb_vec_clear(series, n);

  arb_clear(rest_term);
  arb_clear(tmp);
}

/* Enclose the maximum value of the eigenfunction with by the given
 * coefficients and parameters. The argument vertex indicates which
 * vertex the eigenfunction should be computed from, the boundary
 * which the maximum is computed on is the one opposite to the given
 * vertex.
 *
 * The computations are done by adaptively splitting the interval into
 * smaller intervals and using maximize_series to enclose the maximum
 * on these smaller intervals.
 *
 * TODO: Handle increasing tolerance if required. */
void
maximize(arb_t max, geom_t geom, arb_ptr coefs, slong N, arb_t nu,
         slong vertex, slong prec)
{
  /* Balls are not optimal for storing intervals and we therefore
   * choose to store the intervals as two balls, with radius zero,
   * representing left and right endpoint. We store a list of the
   * intervals we are currently working on and one of the intervals to
   * use for the next iteration. */
  arb_ptr intervals_lower, intervals_upper;
  arb_ptr next_intervals_lower, next_intervals_upper;

  arb_ptr evals;
  arb_t t;
  arf_t max_lower, max_upper, tmp;
  slong n, intervals_len, next_intervals_len;
  slong iterations, max_iterations, splits;

  int done;

  /* Number of terms to use in the Taylor expansion in
   * maximize_series, a good heuristics seems to be to use the same
   * number of terms as in the expansion for the eigenfunction. */
  n = N;

  /* We set the maximum number of iterations to avoid infinite
   * splitting. */
  max_iterations = 30;

  /* The boundary is parameterized with t in the interval [0, 1] and
   * we start with using the whole interval. */
  intervals_len = 1;
  intervals_lower = _arb_vec_init(1);
  intervals_upper = _arb_vec_init(1);
  arb_set_si(intervals_lower + 0, 0);
  arb_set_si(intervals_upper + 0, 1);

  arb_init(t);

  arf_init(max_lower);
  arf_init(max_upper);
  arf_init(tmp);

  done = 0;
  iterations = 0;
  splits = 0;

  while (!done)
  {
    /* Create vector for storing maximum on all the remaining
     * intervals. */
    evals = _arb_vec_init(intervals_len);

    /* Reset values from previous iteration. */
    arb_zero(max);
    arf_zero(max_upper);
    arf_zero(max_lower);

    /* Compute enclosures of the maximum for all intervals using
     * maximize_series on each one. */
    for (slong i = 0; i < intervals_len; i++)
    {
      arb_union(t, intervals_lower + i, intervals_upper + i, prec);

      maximize_series(evals + i, geom, t, coefs, N, nu, n, vertex, prec);

      /* Update lower and upper bounds for the maximum. */
      if (arb_is_finite(evals + i))
      {
        arb_get_abs_ubound_arf(tmp, evals + i, prec);
        arf_max(max_upper, max_upper, tmp);
        arb_get_abs_lbound_arf(tmp, evals + i, prec);
        arf_max(max_lower, max_lower, tmp);
      }
      else
      {
        arf_pos_inf(max_upper);
      }
    }

    /* Set enclosure containing the maximum. */
    arb_set_interval_arf(max, max_lower, max_upper, prec);

    if (arb_rel_error_bits(max) < -5 || iterations > max_iterations)
    {
      /* If the enclosure of the maximum is sufficiently tight or if
       * we have reached the maximum number of iterations we exit. */
      done = 1;
    }
    else
    {
      /* If the enclosure is not tight enough we split all the
       * intervals where the maximum could be obtained. */

      /* Create vertors for storing the intervals for the next
       * iteration, at most we have 2*interval_len ones. */
      next_intervals_lower = _arb_vec_init(2*intervals_len);
      next_intervals_upper = _arb_vec_init(2*intervals_len);
      next_intervals_len = 0;

      for (slong i = 0; i < intervals_len; i++)
      {

        arb_get_abs_ubound_arf(tmp, evals + i, prec);
        if (!(arf_cmp(tmp, max_lower) < 0))
        {
          /* Compute the midpoint of the current interval. */
          arb_add(t, intervals_lower + i, intervals_upper + i, prec);
          arb_div_si(t, t, 2, prec);

          /* To avoid accumulating errors we set the radius of the
           * ball to zero. It is not required that we split exactly at
           * the midpoint so it will not give any problems. */
          mag_zero(arb_radref(t));

          /* Left interval */
          arb_set(next_intervals_lower + next_intervals_len,
                  intervals_lower + i);
          arb_set(next_intervals_upper + next_intervals_len,
                  t);
          /* Right interval */
          arb_set(next_intervals_lower + next_intervals_len + 1,
                  t);
          arb_set(next_intervals_upper + next_intervals_len + 1,
                  intervals_upper + i);

          next_intervals_len += 2;
          splits += 1;
        }
      }

      /* Clear variables for next iteration. */
      _arb_vec_clear(intervals_lower, intervals_len);
      _arb_vec_clear(intervals_upper, intervals_len);
      _arb_vec_clear(evals, intervals_len);

      /* Set variables to use for next iteration. */
      intervals_lower = _arb_vec_init(next_intervals_len);
      intervals_upper = _arb_vec_init(next_intervals_len);
      _arb_vec_set(intervals_lower, next_intervals_lower, next_intervals_len);
      _arb_vec_set(intervals_upper, next_intervals_upper, next_intervals_len);
      _arb_vec_clear(next_intervals_lower, 2*intervals_len);
      _arb_vec_clear(next_intervals_upper, 2*intervals_len);

      intervals_len = next_intervals_len;

      iterations++;
    }
  }

  /* Debug information showing the number of splits and iterations
   * when computing the maximum. */
#ifdef DEBUG
  flint_printf("DEBUG maximize: Splits: %d, Iterations: %d, Maximum: ",
               splits, iterations);
  arb_printn(max, 10, 0);
  flint_printf("\n");
#endif

  _arb_vec_clear(intervals_lower, intervals_len);
  _arb_vec_clear(intervals_upper, intervals_len);
  _arb_vec_clear(evals, intervals_len);

  arb_clear(t);

  arf_clear(max_lower);
  arf_clear(max_upper);
  arf_clear(tmp);
}
