#include "norm.h"

#include "acb_hypgeom.h"
#include "acb_calc.h"

static int
integrand(acb_ptr out, const acb_t inp, void *param_void, slong order,
          slong prec)
{
  acb_ptr nu, mu, coef;
  acb_t tmp;
  slong prec_local;

  acb_init(tmp);

  nu = (acb_ptr)(param_void);
  mu = (acb_ptr)(param_void) + 1;
  coef = (acb_ptr)(param_void) + 2;

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
    acb_hypgeom_legendre_p(tmp, nu, mu, out, 0, prec_local);
    prec_local *= 2;
  } while (!acb_is_finite(tmp) && prec_local <= 16*prec);

  acb_sqr(out, tmp, prec);
  acb_sin(tmp, inp, prec);
  acb_mul(out, out, tmp, prec);
  acb_mul(out, out, coef, prec);

  acb_clear(tmp);

  return 0;
}

/* Compute a lower bound of the norm of the eigenfunction. */
void
integral_norm(arb_t norm, geom_t geom, arb_ptr coefs, int N, arb_t nu,
              slong vertex, slong prec)
{
  acb_t a, b, norm_part;
  acb_ptr params;
  arb_t integral_phi, mu;
  mag_t tol;
  slong prec_local, n;

  params = _acb_vec_init(3);

  acb_init(a);
  acb_init(b);
  acb_init(norm_part);
  arb_init(integral_phi);
  arb_init(mu);

  mag_init(tol);

  arb_set(acb_realref(params + 0), nu); /* params[0] = nu */

  /* The integrals goes from zero to the lower bound for theta.
   * However the function has a branch cut at zero and Arb has problem
   * handling this. We therefore integrate a small distance away from
   * zero. */
  arb_set_str(acb_realref(a), "1e-1", prec);
  arb_set(acb_realref(b), geom->theta_lower[vertex]);

  /* Compute the integral corresponding to phi */
  arb_const_pi(integral_phi, prec);
  arb_div_si(integral_phi, integral_phi, 2, prec);
  geom_get_mu(mu, geom, vertex, 0, prec);
  arb_div(integral_phi, integral_phi, mu, prec);
  arb_neg(integral_phi, integral_phi);

  /* The norms of the basis functions are rapidly decreasing so it is
     enough to only compute the first few to get a good lower
     bound. */

  if (N < 4)
  {
    n = N;
  }
  else
  {
    n = 4;
  }

  for (slong i = 0; i < n; i++)
  {
    /* params[1] = mu */
    geom_get_mu(acb_realref(params + 1), geom, vertex, i, prec);
    /* params[2] = sqr(coefs[i]) */
    arb_sqr(acb_realref(params + 2), coefs + i, prec);

    /* Compute the integral corresponding to theta */
    prec_local = prec;
    do
    {
      acb_calc_integrate(norm_part, integrand, params, a, b, prec, tol, NULL,
                         prec_local);
      prec_local *= 2;
    } while (!acb_is_finite(norm_part));

    /* Add to the norm */
    arb_add(norm, norm, acb_realref(norm_part), prec);
  }
  /* Multiply with the integral corresponding to phi */
  arb_mul(norm, norm, integral_phi, prec);

  arb_sqrt(norm, norm, prec);

  _acb_vec_clear(params, 3);

  acb_clear(a);
  acb_clear(b);
  acb_clear(norm_part);
  arb_clear(integral_phi);
  arb_clear(mu);

  mag_clear(tol);
}
