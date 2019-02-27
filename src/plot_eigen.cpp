#include "plot_eigen.h"

#include "generate-matrix.h"
#include "eigenfunction.h"
#include "arb_poly.h"

void
plot_eigen(geom_t geom, arb_ptr coefs, slong N, arb_t nu, slong n, slong minces,
           slong vertex, slong prec) {
  arb_ptr z, phi, poly;
  arb_t t, mid, theta, res, rest_term;
  arf_t a, b;

  vertex = 0;
  minces = 250;

  arb_init(t);
  arb_init(mid);
  arb_init(theta);
  arb_init(res);
  arb_init(rest_term);

  arf_init(a);
  arf_init(b);

  if (n < 0)
  {
    z = _arb_vec_init(1);
    phi = _arb_vec_init(1);
    poly = _arb_vec_init(1);

    /* Non rigorous plot, only plot on points */
    for (slong i = 0; i < minces; i++)
    {
      arf_set_si(a, i);
      arf_div_si(a, a, minces, prec, ARF_RND_DOWN);
      arf_set_si(b, i + 1);
      arf_div_si(b, b, minces, prec, ARF_RND_UP);

      arb_set_interval_arf(t, a, b, prec);
      arb_set_arf(mid, arb_midref(t));

      /* Compute the first coefficient of the Taylor polynomial of the
       * eigenfunction at the midpoint of t */
      parametrization(z, phi, geom, mid, 1, vertex, prec);
      eigenfunction_series(poly, geom, coefs, N, z, phi, nu, vertex, 1, prec);
      arb_set(res, poly + 0);

      /* To keep the same output format as the rigorous case we print
       * both the left and right endpoint of the midpoint, even though
       * they are the same. */
      arb_get_interval_arf(a, b, mid, prec);
      arf_printd(a, 10);
      flint_printf(" ");
      arf_printd(b, 10);
      flint_printf(" ");

      arb_get_interval_arf(a, b, res, prec);
      arf_printd(a, 10);
      flint_printf(" ");
      arf_printd(b, 10);
      flint_printf("\n");
    }

    _arb_vec_clear(z, 1);
    _arb_vec_clear(phi, 1);
    _arb_vec_clear(poly, 1);
  }
  else
  {
    z = _arb_vec_init(n + 1);
    phi = _arb_vec_init(n + 1);
    poly = _arb_vec_init(n + 1);

    for (slong i = 0; i < minces; i++)
    {
      arf_set_si(a, i);
      arf_div_si(a, a, minces, prec, ARF_RND_DOWN);
      arf_set_si(b, i + 1);
      arf_div_si(b, b, minces, prec, ARF_RND_UP);

      arb_set_interval_arf(t, a, b, prec);

      if (n > 0)
      {
        /* Compute the Taylor polynomial of the eigenfunction at the midpoint of t */
        arb_set_arf(mid, arb_midref(t));
        parametrization(z, phi, geom, mid, n, vertex, prec);
        eigenfunction_series(poly, geom, coefs, N, z, phi, nu, vertex, n, prec);
      }

      /* Enclose the Taylor polynomial evaluated at t - mid(t) */
      arb_sub(mid, t, mid, prec);
      _arb_poly_evaluate(res, poly, n, mid, prec);

      /* Compute the rest term of the Taylor expansion */
      parametrization(z, phi, geom, t, n + 1, vertex, prec);
      eigenfunction_series(poly, geom, coefs, N, z, phi, nu, vertex, n + 1, prec);

      arb_pow_ui(rest_term, mid, n, prec);
      arb_mul(rest_term, rest_term, poly + n, prec);

      /* Add the rest term to the maximum */
      arb_add(res, res, rest_term, prec);

      arb_get_interval_arf(a, b, t, prec);
      arf_printd(a, 10);
      flint_printf(" ");
      arf_printd(b, 10);
      flint_printf(" ");

      arb_get_interval_arf(a, b, res, prec);
      arf_printd(a, 10);
      flint_printf(" ");
      arf_printd(b, 10);
      flint_printf("\n");
    }

    _arb_vec_clear(z, n + 1);
    _arb_vec_clear(phi, n + 1);
    _arb_vec_clear(poly, n + 1);

  }

  arb_clear(t);
  arb_clear(mid);
  arb_clear(theta);
  arb_clear(res);
  arb_clear(rest_term);

  arf_clear(a);
  arf_clear(b);
}
