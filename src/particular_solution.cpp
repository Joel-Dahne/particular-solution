#include "particular_solution.h"

#include "sigma.h"
#include "enclose.h"
#include "plot_eigen.h"

void
particular_solution_opt_init(particular_solution_opt_t options)
{
  arb_init(options->prec_factor);
  arb_init(options->tol_relative);
}

void
particular_solution_opt_default(particular_solution_opt_t options)
{
  arb_set_d(options->prec_factor, 1.2);
  arb_set_d(options->tol_relative, 1e-5);
  options->N_beg = 4;
  options->N_end = 16;
  options->N_step = 2;
  options->output = 0;
}

void
particular_solution_opt_clear(particular_solution_opt_t options)
{
  arb_clear(options->prec_factor);
  arb_clear(options->tol_relative);
}

void
particular_solution_enclosure(arb_t nu_enclosure, geom_t geometry,
                              particular_solution_opt_t options, slong prec)
{
  arb_ptr coefs;
  arb_t nu, tol, tmp;
  points_t points;

  arb_init(nu);
  arb_init(tol);
  arb_init(tmp);

  for (slong N = options->N_beg; N <= options->N_end; N += options->N_step)
  {
    /* Compute tolerance and precision to use */
    arb_get_rad_arb(tol, nu_enclosure);
    arb_mul(tol, tol, options->tol_relative, prec);

    /* FIXME: Use the size of nu in the computations for the precision */
    arb_log_base_ui(tmp, tol, 2, prec);
    arb_mul(tmp, tmp, options->prec_factor, prec);
    arb_neg(tmp, tmp);

    if (arf_get_si(arb_midref(tmp), ARF_RND_CEIL) > prec)
    {
      prec = arf_get_si(arb_midref(tmp), ARF_RND_CEIL);
    }

    mpfr_set_default_prec(prec);

    /* Recompute variables to new precision */
    geom_compute(geometry, prec);

    /* Initiate new variables */
    coefs = _arb_vec_init(N);

    points_init(points, 2*N, 2*N);

    boundary(points, geometry, prec);
    interior(points, geometry, prec);

    /* Find the value of nu that minimizes sigma */
    minimize_sigma(nu, geometry, points, N, nu_enclosure, tol, prec);

    /* Find the coefficients of the expansion */
    coefs_sigma(coefs, geometry, points, N, nu, prec);

    /* Compute an enclosure of the eigenvalue. To be sure to get correct
       output of the eigenvalue the output of it is handled inside the
       function enclose. */
    enclose(nu_enclosure, geometry, coefs, N, nu, prec);

    /* Print information */
    if (options->output != 0)
    {
      flint_printf("%i ", N);

      /* Compute bounds for the eigenvalue */
      arb_add_si(tmp, nu_enclosure, 1, prec);
      arb_mul(tmp, tmp, nu_enclosure, prec);

      if (options->output == 1)
      {
        arb_printn(tmp, (slong)ceil(prec*log10(2)), 0);
        flint_printf("\n");
      }
      else if (options->output == 2)
      {
        arb_printn(nu_enclosure, (slong)ceil(prec*log10(2)), 0);
        flint_printf("\n");
      }
      else if (options->output == 3)
      {
        arf_printd(arb_midref(tmp), (slong)ceil(prec*log10(2)));
        flint_printf("\n");
      }
      else if (options->output == 4)
      {
        arf_printd(arb_midref(nu_enclosure), (slong)ceil(prec*log10(2)));
        flint_printf("\n");
      }
      else if (options->output == 5)
      {
        flint_printf(" %e\n", mag_get_d(arb_radref(tmp)));
      }
      else if (options->output == 6)
      {
        flint_printf(" %e\n", mag_get_d(arb_radref(nu_enclosure)));
      }
      else if (options->output == 7)
      {
        flint_printf("\n");
        for (slong i = 0; i < N; i++)
        {
          arf_printd(arb_midref(coefs + i), (slong)ceil(prec*log10(2)));
          flint_printf("\n");
        }
      }
      else if (options->output == 8)
      {
        arf_printd(arb_midref(nu), (slong)ceil(prec*log10(2)));
        flint_printf(" %i\n", 500);
        plot_eigen(geometry, coefs, N, nu, 500, prec);
      }


    }

    _arb_vec_clear(coefs, N);

    points_clear(points);
  }

  arb_clear(nu);
  arb_clear(tol);
  arb_clear(tmp);
}
