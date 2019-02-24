#include "eigenfunction.h"

#include "arb_hypgeom.h"

/* Compute the Taylor expansion of the Legendre P function evaluated
 * at the Taylor expansion given by z with parameters nu and mu. The
 * order of the computed expansion is n - 1, i.e. the number of terms
 * is equal to n.
 *
 * For the computations it uses a fourth order recusion formula for
 * the terms in the expansion which can be deduced from the
 * differential equations defining the Legendre functions. The first
 * four terms are computed explicitly. */
static void
legendre_p_series(arb_ptr res, arb_ptr z, arb_t nu, arb_t mu, slong n,
                  slong prec)
{
  arb_ptr legendre, tmp_poly;
  arb_t z2, nu2, mu2, onemz2, tmp1, tmp2;
  slong prec_local, k;

  legendre = _arb_vec_init(n);
  tmp_poly = _arb_vec_init(n);

  arb_init(z2);
  arb_init(nu2);
  arb_init(mu2);
  arb_init(onemz2);
  arb_init(tmp1);
  arb_init(tmp2);

  /* z2 = z^2 */
  arb_sqr(z2, z, prec);
  /* nu2 = nu^2 */
  arb_sqr(nu2, nu, prec);
  /* mu2 = mu^2 */
  arb_sqr(mu2, mu, prec);
  /* onemz2 = 1 - z^2 */
  arb_neg(onemz2, z2);
  arb_add_si(onemz2, onemz2, 1, prec);

  if (n > 0)
  {
    /* a_0 */
    prec_local = prec;
    do {
      arb_hypgeom_legendre_p(legendre + 0, nu, mu, z, 0, prec_local);
      prec_local = 2*prec_local;
    } while (!arb_is_finite(legendre + 0) && prec_local < 16*prec);

    if (n > 1)
    {
      /* a_1 */
      arb_add_si(tmp1, nu, 1, prec);
      prec_local = prec;
      do {
        arb_hypgeom_legendre_p(legendre + 1, tmp1, mu, z, 0, prec_local);
        prec_local = 2*prec_local;
      } while (!arb_is_finite(legendre + 1) && prec_local < 16*prec);

      arb_sub(tmp1, tmp1, mu, prec);
      arb_mul(legendre + 1, legendre + 1, tmp1, prec);

      arb_add_si(tmp1, nu, 1, prec);
      arb_mul(tmp1, tmp1, z, prec);
      arb_mul(tmp1, tmp1, legendre + 0, prec);
      arb_sub(legendre + 1, tmp1, legendre + 1, prec);

      arb_div(legendre + 1, legendre + 1, onemz2, prec);

      if (n > 2)
      {
        /* a_2 */
        arb_mul_si(tmp1, z, 2, prec);
        arb_mul(legendre + 2, tmp1, legendre + 1, prec);

        arb_div(tmp1, mu2, onemz2, prec);
        arb_add_si(tmp2, nu, 1, prec);
        arb_mul(tmp2, tmp2, nu, prec);
        arb_sub(tmp1, tmp2, tmp1, prec);
        arb_submul(legendre + 2, tmp1, legendre + 0, prec);

        arb_div(legendre + 2, legendre + 2, onemz2, prec);
        arb_div_si(legendre + 2, legendre + 2, 2, prec);

        if (n > 3)
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
          for (slong i = 4; i < n; i++)
          {
            k = i - 4;

            /* Term for i */
            arb_add_si(tmp1, nu, k + 1, prec);
            arb_add_si(tmp2, nu, -k, prec);
            arb_mul(tmp1, tmp1, tmp2, prec);
            arb_mul(legendre + k + 4, tmp1, legendre + k, prec);

            /* Term for i + 1 */
            arb_add(tmp1, nu2, nu, prec);
            arb_neg(tmp1, tmp1);
            arb_add_si(tmp1, tmp1, 5*k + 3, prec);
            arb_div_si(tmp1, tmp1, 2, prec);
            arb_add_si(tmp1, tmp1, k*k, prec);
            arb_mul_si(tmp1, tmp1, -4, prec);
            arb_mul(tmp1, tmp1, z, prec);
            arb_addmul(legendre + k + 4, tmp1, legendre + k + 1, prec);

            /* Term for i + 2 */
            arb_add(tmp1, nu, nu2, prec);
            arb_add_si(tmp1, tmp1, -6*k*k - 24*k - 24, prec);
            arb_mul(tmp1, tmp1, z2, prec);
            arb_add(tmp1, tmp1, mu2, prec);
            arb_sub(tmp1, tmp1, nu2, prec);
            arb_sub(tmp1, tmp1, nu, prec);
            arb_add_si(tmp1, tmp1, 2*k*k + 8*k + 8, prec);
            arb_addmul(legendre + k + 4, tmp1, legendre + k + 2, prec);

            /* Term for i + 3 */
            arb_add_si(tmp1, z, 1, prec);
            arb_mul(tmp1, tmp1, z, prec);
            arb_add_si(tmp2, z, -1, prec);
            arb_mul(tmp1, tmp1, tmp2, prec);
            arb_mul_si(tmp1, tmp1, -2*(k + 3)*(2*k + 5), prec);
            arb_addmul(legendre + k + 4, tmp1, legendre + k + 3, prec);

            /* Factor */
            arb_add_si(tmp1, z, -1, prec);
            arb_sqr(tmp1, tmp1, prec);
            arb_add_si(tmp2, z, 1, prec);
            arb_sqr(tmp2, tmp2, prec);
            arb_mul(tmp1, tmp1, tmp2, prec);
            arb_mul_si(tmp1, tmp1, (k + 4)*(k + 3), prec);
            arb_div(legendre + k + 4, legendre + k + 4, tmp1, prec);
          }
        }
      }
    }
  }

  /* Compose the Taylor series for the Legendre function with that of
   * z */
  _arb_vec_set(tmp_poly, z, n);
  arb_zero(tmp_poly);
  _arb_poly_compose_series(res, legendre, n, tmp_poly, n, n, prec);

  _arb_vec_clear(legendre, n);
  _arb_vec_clear(tmp_poly, n);

  arb_clear(z2);
  arb_clear(nu2);
  arb_clear(mu2);
  arb_clear(onemz2);
  arb_clear(tmp1);
  arb_clear(tmp2);
}

/* Compute the value of the basis functions used in the expansion for
 * the eigenfunction at the point given by theta and phi with the
 * paremeters nu and mu. */
void
eigenfunction_basis(arb_t res, arb_t theta, arb_t phi, arb_t nu, arb_t mu,
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

/* Compute the Taylor expansion of the basis functions used in the
 * expansion for the eigenfunction using the Taylor expansions for z
 * and phi as input with the parameters nu and mu. The order of the
 * computed expansion is n - 1, i.e. the number of terms is equal to
 * n. */
void
eigenfunction_basis_series(arb_ptr res, arb_ptr z, arb_ptr phi, arb_t nu,
                           arb_ptr mu, slong n, slong prec)
{
  arb_ptr tmp1, tmp2;

  tmp1 = _arb_vec_init(n);
  tmp2 = _arb_vec_init(n);

  _arb_vec_scalar_mul(tmp1, phi, n, mu, prec);
  _arb_poly_sin_series(tmp2, tmp1, n, n, prec);

  legendre_p_series(tmp1, z, nu, mu, n, prec);

  _arb_poly_mullow(res, tmp1, n, tmp2, n, n, prec);

  _arb_vec_clear(tmp1, n);
  _arb_vec_clear(tmp2, n);
}

/* Compute the value of the eigenfunction function given by a sum of N
 * basis functions with the given coefficients at the point given by
 * theta and phi using the parameter nu. The argument vertex indicates
 * which vertex the eigenfunction should be computed from, this is
 * used for deciding the value of the parameter mu for the basis
 * functions. */
void
eigenfunction(arb_t res, geom_t geom, arb_ptr coefs, slong N, arb_t theta,
              arb_t phi, arb_t nu, slong vertex, slong prec)
{
  arb_t mu, term;

  arb_init(mu);
  arb_init(term);

  arb_zero(res);

  for (slong i = 0; i < N; i++)
  {
    geom_get_mu(mu, geom, vertex, i, prec);

    eigenfunction_basis(term, theta, phi, nu, mu, prec);

    arb_addmul(res, term, coefs + i, prec);
  }

  arb_clear(mu);
  arb_clear(term);
}

/* Compute the Taylor expansion of the eigenfunction function given by
 * a sum of N basis functions with the given coefficients at the point
 * given by z and phi using the parameter nu. The argument vertex
 * indicates which vertex the eigenfunction should be computed from,
 * this is used for deciding the value of the parameter mu for the
 * basis functions. The order of the computed expansion is n - 1, i.e.
 * the number of terms is equal to n. */
void
eigenfunction_series(arb_ptr res, geom_t geom, arb_ptr coefs, slong N,
                     arb_ptr z, arb_ptr phi, arb_t nu, slong vertex, slong n,
                     slong prec)
{
  arb_ptr term;
  arb_t mu;

  term = _arb_vec_init(n);

  arb_init(mu);

  _arb_vec_zero(res, n);

  for (slong i = 0; i < N; i++)
  {
    geom_get_mu(mu, geom, vertex, i, prec);

    eigenfunction_basis_series(term, z, phi, nu, mu, n, prec);

    _arb_vec_scalar_addmul(res, term, n, coefs + i, prec);
  }

  _arb_vec_clear(term, n);

  arb_clear(mu);
}
