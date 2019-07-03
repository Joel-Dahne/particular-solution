#include "particular_solution.h"

#include "sigma.h"
#include "enclose.h"
#include "plot_eigen.h"
#include "time.h"
#include <math.h>

void
particular_solution_enclosure(arb_t nu_enclosure, geom_t geometry,
                              options_t options, slong prec)
{
  arb_ptr* coefs;
  arb_t nu, tol, tmp;
  points_t points;
  clock_t start_time, sigma_time, enclosure_time;
  slong vertex;

  arb_init(nu);
  arb_init(tol);
  arb_init(tmp);

  coefs = (arb_ptr*)calloc(3, sizeof(arb_ptr));

  for (slong N = options->N_beg; N <= options->N_end; N += options->N_step)
  {
    /* Compute tolerance and precision to use */
    arb_get_rad_arb(tol, nu_enclosure);
    arb_set_d(tmp, options->tol_relative);
    arb_mul(tol, tol, tmp, prec);

    /* FIXME: Use the size of nu in the computations for the precision */
    arb_log_base_ui(tmp, tol, 2, prec);
    arb_neg(tmp, tmp);

    if (arf_get_si(arb_midref(tmp), ARF_RND_CEIL)*options->prec_factor > prec)
    {
      prec = arf_get_si(arb_midref(tmp), ARF_RND_CEIL)*options->prec_factor;
    }

    /* Recompute variables to take into account the new precision.
     * Certain parts of the program works with higher intermediary
     * precision, to not have these parts be limited by the precision
     * of these pre computed variables we compute them to a higher
     * precision than we are currently working with. Testing has
     * showed that using double the current precision seems to work
     * well. */
    geom_compute(geometry, 2*prec);

    /* Initiate new variables */
    for (slong i = 0; i < 3; i++)
    {
      coefs[i] = _arb_vec_init(N);
    }

    points_init(points, geometry, 2*N, 2*N);

    boundary(points, geometry, prec);
    interior(points, geometry, prec);

    start_time = clock();

    /* Find the value of nu that minimizes sigma */
    minimize_sigma(nu, geometry, points, N, nu_enclosure, tol, prec);

    /* Find the coefficients of the expansion */
    coefs_sigma(coefs, geometry, points, N, nu, prec);
    sigma_time = clock();

    /* Compute an enclosure of the eigenvalue. */
    /* FIXME: At the moment enclose only supports using a single
     * expansion from one vertex and doesn't handle expansions from
     * several vertices. We choose which vertex to use by taking the
     * first vertex we have an expansion from. This should be changed
     * once enclose supports expansions from several vertices. */
    if (geometry->vertices[0])
    {
      vertex = 0;
    }
    else if (geometry->vertices[1])
    {
      vertex = 1;
    }
    else
    {
      vertex = 2;
    }

    enclose(nu_enclosure, geometry, coefs, N, nu, vertex, prec);
    enclosure_time = clock();

    /* Print information */
    if (options->output != 0
        && (!options->output_final || (options->N_end == N)))
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
          arf_printd(arb_midref(coefs[0] + i), (slong)ceil(prec*log10(2)));
          flint_printf("\n");
        }
      }
      else if (options->output == 8)
      {
        arf_printd(arb_midref(nu), (slong)ceil(prec*log10(2)));
        flint_printf(" %i\n", 500);
        plot_eigen(geometry, coefs[0], N, nu, options->plot_n, 250, 0, prec);
      }
      else if (options->output == 9)
      {
        arf_printd(arb_midref(nu), (slong)ceil(prec*log10(2)));
        flint_printf("\n");
      }

    }

    if (options->output_time
         && (!options->output_final || (options->N_end == N)))
    {
      flint_printf("Time computing minimum: %f\n",
                   ((double)(sigma_time - start_time))/CLOCKS_PER_SEC);
      flint_printf("Time computing enclosure: %f\n",
                   ((double)(enclosure_time - sigma_time))/CLOCKS_PER_SEC);
    }

    for (slong i = 0; i < 3; i++)
    {
      _arb_vec_clear(coefs[i], N);
    }

    points_clear(points);
  }

  free(coefs);
  coefs = NULL;

  arb_clear(nu);
  arb_clear(tol);
  arb_clear(tmp);
}
