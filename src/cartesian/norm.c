#include "norm.h"

#include "acb_hypgeom.h"
#include "acb_calc.h"

static int
integrand(acb_ptr out, const acb_t inp, void *param_void, slong order,
          slong prec)
{
  acb_ptr sqrtlambda, nu;
  acb_t tmp;
  slong prec_local;

  acb_init(tmp);

  sqrtlambda = (acb_ptr)(param_void);
  nu = (acb_ptr)(param_void) + 1;

  if (order > 1)
  {
    flint_abort();
  }

  acb_mul(out, sqrtlambda, inp, prec);

  if (order == 1)
  {
    if (arb_contains_negative(acb_realref(out)))
    {
      acb_indeterminate(out);
      return 0;
    }
  }

  prec_local = prec;
  do
  {
    acb_hypgeom_bessel_j(tmp, nu, out, prec);
    prec_local *= 2;
  } while (!acb_is_finite(tmp) && prec_local <= 16*prec);

  acb_sqr(out, tmp, prec);
  acb_mul(out, out, inp, prec);

  acb_clear(tmp);

  return 0;
}

void
integral_norm(arb_t norm, geom_t geom, arb_ptr coefs, int N, arb_t lambda,
              slong vertex, slong prec)
{
  acb_t a, b, norm_part;
  acb_ptr params;
  arb_t integral_theta, nu, tmp;
  mag_t tol;
  slong prec_local, n;

  params = _acb_vec_init(2);

  acb_init(a);
  acb_init(b);
  acb_init(norm_part);
  arb_init(integral_theta);
  arb_init(nu);
  arb_init(tmp);

  mag_init(tol);

  arb_sqrt(acb_realref(params + 0), lambda, prec); /* params[0] = sqrt(lambda) */

  /* The integrals goes from zero to the lower bound for theta.
   * However the function has a branch cut at zero and Arb has problem
   * handling this. We therefore integrate a small distance away from
   * zero. */
  arb_set_d(acb_realref(a), 0.5);
  arb_set(acb_realref(b), geom->r_lower);

  /* Compute the integral corresponding to theta */
  arb_const_pi(integral_theta, prec);
  arb_div_si(integral_theta, integral_theta, 2, prec);
  geom_get_mu(nu, geom, vertex, 0, prec);
  arb_div(integral_theta, integral_theta, nu, prec);

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
    /* params[1] = nu */
    geom_get_mu(acb_realref(params + 1), geom, vertex, i, prec);

    /* Compute the integral corresponding to theta */
    prec_local = prec;
    do
    {
      acb_calc_integrate(norm_part, integrand, params, a, b, prec, tol, NULL,
                         prec_local);

      prec_local *= 2;
    } while (!acb_is_finite(norm_part));

    /* Multiply with the coefficient squared and add to the integral */
    arb_sqr(tmp, coefs + i, prec);
    arb_addmul(norm, acb_realref(norm_part), tmp, prec);
  }
  /* Multiply with the integral corresponding to phi */
  arb_mul(norm, norm, integral_theta, prec);

  arb_sqrt(norm, norm, prec);

#ifdef DEBUG
  flint_printf("DEBUG integral_norm: ");
  arb_printn(norm, 10, 0);

  flint_printf("\n");
#endif

  _acb_vec_clear(params, 2);

  acb_clear(a);
  acb_clear(b);
  acb_clear(norm_part);
  arb_clear(integral_theta);
  arb_clear(nu);
  arb_clear(tmp);

  mag_clear(tol);
}
