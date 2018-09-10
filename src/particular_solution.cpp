#include "particular_solution.h"

#include "sigma.h"
#include "enclose.h"

int
index_function_odd(int k) {
  return 2*k + 1;
}

int
index_function_all(int k) {
  return k + 1;
}

void
particular_solution_opt_init(particular_solution_opt_t options)
{
  options->prec_factor = 1.2;
  options->N_beg = 4;
  options->N_end = 16;
  options->N_step = 2;
  options->verbose = 0;
  options->index_function = index_function_all;
}

void
particular_solution_enclosure(geom_t geometry, int angles_coefs[],
                              mpfr_t mu0, mpfr_t nu_low, mpfr_t nu_upp, mpfr_t tol,
                              particular_solution_opt_t options)
{
  mpfr_t *coefs;
  mpfr_t nu;
  points_t points;

  mpfr_init(nu);

  for (int N = options->N_beg; N <= options->N_end; N += options->N_step)
  {
    /* Compute precision to use */

    /* Initiate and update variables to the precision */
    coefs = new mpfr_t[N];
    for (int i = 0; i < N; i++) {
      mpfr_init(coefs[i]);
    }

    points_init(points, 2*N, 2*N);

    boundary(points, geometry);
    interior(points, geometry);

    /* Find the value of nu that minimizes sigma */
    minimize_sigma(nu, points, N, nu_low, nu_upp, mu0, tol, options->index_function);

    mpfr_printf("%i ", N);
    fflush(stdout);

    /* Find the coefficients of the expansion */
    coefs_sigma(coefs, points, N, nu, mu0, options->index_function);

    /* Compute an enclosure of the eigenvalue. To be sure to get correct
       output of the eigenvalue the output of it is handled inside the
       function enclose. */
    enclose(nu_low, nu_upp, angles_coefs, coefs, N, nu, options->index_function, 4);
    mpfr_printf("\n");

    for (int i = 0; i < N; i++) {
      mpfr_clear(coefs[i]);
    }
    delete [] coefs;

    points_clear(points);
  }

  mpfr_clear(nu);
}
