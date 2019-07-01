#include "plot_sigma.h"

#include "sigma.h"

void
plot_sigma(arf_t inf, arf_t sup, geom_t geometry, slong num_points,
           particular_solution_opt_t options, slong prec) {
  arb_t nu, res, step;
  points_t points;

  mpfr_set_default_prec(prec);

  arb_init(nu);
  arb_init(res);
  arb_init(step);

  geom_compute(geometry, prec);

  arb_set_arf(step, sup);
  arb_sub_arf(step, step, inf, prec);
  arb_div_si(step, step, num_points, prec);

  for (slong N = options->N_beg; N <= options->N_end; N += options->N_step)
  {
    points_init(points, geometry, 2*N, 2*N);

    boundary(points, geometry, prec);
    interior(points, geometry, prec);

    flint_printf("%i %i\n", N, num_points);

    for (slong i = 0; i < num_points; i++)
    {
      arb_set_arf(nu, inf);
      arb_addmul_si(nu, step, i, prec);

      sigma(res, geometry, points, N, nu, prec);

      arf_printd(arb_midref(nu), (slong)ceil(prec*log10(2)));
      arf_printd(arb_midref(res), (slong)ceil(prec*log10(2)));
      flint_printf("\n");
    }

    points_clear(points);
  }

  arb_clear(nu);
  arb_clear(res);
  arb_clear(step);
}
