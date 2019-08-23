#include "eigenfunction.h"

#include "arb_hypgeom.h"

/* Compute the Taylor expansion of the Bessel J function evaluated at
 * the Taylor expansion given by x with parameter nu. The order of the
 * computed expansion is n - 1, i.e. the number of terms is equal to
 * n.
 *
 * For the computations it uses ... */
static void
bessel_j_series(arb_ptr res, arb_t nu, arb_ptr x, slong len, slong n,
                slong prec)
{
  arb_ptr bessel, tmp_poly;
  arb_t nu_tmp, bessel_tmp, tmp;
  slong k;

  bessel = _arb_vec_init(n);
  tmp_poly = _arb_vec_init(n);

  arb_init(nu_tmp);
  arb_init(bessel_tmp);
  arb_init(tmp);

  if (n > 0)
  {
    /* a_0 */
    arb_hypgeom_bessel_j(bessel + 0, nu, x, prec);

    if (n > 1)
    {
      /* a_1 */
      arb_sub_si(nu_tmp, nu, 1, prec);
      arb_hypgeom_bessel_j(bessel + 1, nu_tmp, x, prec);

      arb_add_si(nu_tmp, nu, 1, prec);
      arb_hypgeom_bessel_j(bessel_tmp, nu_tmp, x, prec);

      arb_sub(bessel + 1, bessel + 1, bessel_tmp, prec);
      arb_div_si(bessel + 1, bessel + 1, 2, prec);

      if (n > 2)
      {
        /* a_2 */
        arb_mul_si(bessel + 2, bessel + 0, -2, prec);

        arb_sub_si(nu_tmp, nu, 2, prec);
        arb_hypgeom_bessel_j(bessel_tmp, nu_tmp, x, prec);
        arb_add(bessel + 2, bessel + 2, bessel_tmp, prec);

        arb_add_si(nu_tmp, nu, 2, prec);
        arb_hypgeom_bessel_j(bessel_tmp, nu_tmp, x, prec);
        arb_add(bessel + 2, bessel + 2, bessel_tmp, prec);

        arb_div_si(bessel + 2, bessel + 2, 8, prec);

        if (n > 3)
        {
          /* a_3 */
          arb_mul_si(bessel + 3, bessel + 1, -6, prec);

          arb_sub_si(nu_tmp, nu, 3, prec);
          arb_hypgeom_bessel_j(bessel_tmp, nu_tmp, x, prec);
          arb_add(bessel + 3, bessel + 3, bessel_tmp, prec);

          arb_add_si(nu_tmp, nu, 3, prec);
          arb_hypgeom_bessel_j(bessel_tmp, nu_tmp, x, prec);
          arb_sub(bessel + 3, bessel + 3, bessel_tmp, prec);

          arb_div_si(bessel + 3, bessel + 3, 48, prec);

          /* Compute the remaining terms using the recurrence relation */
          for (slong i = 4; i < n; i++)
          {
            k = i - 4;

            /* a_k term*/
            arb_set(bessel + k + 4, bessel + k);

            /* a_{k + 1} term*/
            arb_mul_si(tmp, x, 2, prec);

            arb_addmul(bessel + k + 4, bessel + k + 1, tmp, prec);

            /* a_{k + 2} term*/
            arb_sqr(tmp, x, prec);
            arb_submul(tmp, nu, nu, prec);
            arb_add_si(tmp, tmp, k*k + 4*k + 4, prec);

            arb_addmul(bessel + k + 4, bessel + k + 2, tmp, prec);

            /* a_{k + 3} term*/
            arb_mul_si(tmp, x, 2*k*k + 11*k + 15, prec);

            arb_addmul(bessel + k + 4, bessel + k + 3, tmp, prec);

            /* Factor */
            arb_sqr(tmp, x, prec);
            arb_mul_si(tmp, tmp, k*k + 7*k + 12, prec);

            arb_div(bessel + k + 4, bessel + k + 4, tmp, prec);
            arb_neg(bessel + k + 4, bessel + k + 4);
          }
        }
      }
    }
  }

  /* Compose the Taylor series for the Legendre function with that of
   * z */
  _arb_vec_set(tmp_poly, x, n);
  arb_zero(tmp_poly);
  _arb_poly_compose_series(res, bessel, n, tmp_poly, n, n, prec);

  _arb_vec_clear(bessel, n);
  _arb_vec_clear(tmp_poly, n);

  arb_clear(nu_tmp);
  arb_clear(bessel_tmp);
  arb_clear(tmp);
}

/* Compute the value of the basis functions used in the expansion for
 * the eigenfunction at the point given by angle theta and radii r
 * with the parameters lambda an nu. */
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

/* Compute the Taylor expansion of the basis functions used in the
 * expansion for the eigenfunction at the point given by angle theta
 * and radii r with the parameters lambda an nu. The order of the
 * computed expansion is n - 1, i.e. the number of terms is equal to
 * n. */
void
eigenfunction_basis_series(arb_ptr res, arb_ptr r, arb_ptr theta, arb_t lambda,
                           arb_ptr nu, slong n, slong prec)
{
  arb_ptr tmp1, tmp2;
  arb_t sqrtlambda;

  tmp1 = _arb_vec_init(n);
  tmp2 = _arb_vec_init(n);

  arb_init(sqrtlambda);

  arb_sqrt(sqrtlambda, lambda, prec);
  _arb_vec_scalar_mul(tmp1, r, n, sqrtlambda, prec);
  bessel_j_series(tmp2, nu, tmp1, n, n, prec);

  _arb_vec_scalar_mul(tmp1, theta, n, nu, prec);
  _arb_poly_sin_series(tmp1, tmp1, n, n, prec);

  _arb_poly_mullow(res, tmp1, n, tmp2, n, n, prec);

  arb_clear(sqrtlambda);

  /////
  //bessel_j_series(res, nu, r, n, n, prec);

  _arb_vec_clear(tmp1, n);
  _arb_vec_clear(tmp2, n);
}

/* Compute the value of the eigenfunction function given by a sum of N
 * basis functions with the given coefficients at the point given by
 * theta and r using the parameter lambda. */
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

/* Compute the Taylor expansion of the eigenfunction function given by
 * a sum of N basis functions with the given coefficients at the point
 * given by theta and r using the parameter lambda. The order of the
 * computed expansion is n - 1, i.e. the number of terms is equal to
 * n. */
void
eigenfunction_series(arb_ptr res, geom_t geom, arb_ptr coefs, slong N,
                     arb_ptr r, arb_ptr theta, arb_t lambda, slong vertex,
                     slong n, slong prec)
{
  arb_ptr term;
  arb_t nu;

  term = _arb_vec_init(n);

  arb_init(nu);

  _arb_vec_zero(res, n);

  for (slong i = 0; i < N; i++)
  {
    geom_get_mu(nu, geom, vertex, i, prec);

    eigenfunction_basis_series(term, r, theta, lambda, nu, n, prec);

    _arb_vec_scalar_addmul(res, term, n, coefs + i, prec);
  }

  _arb_vec_clear(term, n);

  arb_clear(nu);
}
