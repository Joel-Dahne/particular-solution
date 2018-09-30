#include "plot_sigma.h"

#include "sigma.h"

void
plot_sigma(arf_t inf, arf_t sup, geom_t geometry,
           particular_solution_opt_t options, slong prec) {
  arb_t nu, mu0, step;
  mpfr_t res;
  points_t points;
  slong num_points;

  mpfr_set_default_prec(prec);

  arb_init(nu);
  arb_init(mu0);
  arb_init(step);

  mpfr_init(res);

  num_points = 500;

  geom_compute(geometry, prec);

  arb_set_fmpq(mu0, geometry->angles, prec);
  arb_inv(mu0, mu0, prec);
  arb_neg(mu0, mu0);

  arb_set_arf(step, sup);
  arb_sub_arf(step, step, inf, prec);
  arb_div_si(step, step, num_points, prec);

  for (slong N = options->N_beg; N <= options->N_end; N += options->N_step)
  {
    points_init(points, 2*N, 2*N);

    boundary(points, geometry, prec);
    interior(points, geometry, prec);

    flint_printf("%i %i\n", N, num_points);

    for (slong i = 0; i < num_points; i++)
    {
      arb_set_arf(nu, inf);
      arb_addmul_si(nu, step, i, prec);

      sigma(res, points, N, nu, mu0, options->index_function, prec);

      arf_printd(arb_midref(nu), (slong)ceil(prec*log10(2)));
      mpfr_printf(" %Re\n", res);
    }

    points_clear(points);
  }

  arb_clear(nu);
  arb_clear(mu0);
  arb_clear(step);

  mpfr_clear(res);
}
