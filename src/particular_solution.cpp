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
  options->tol_relative = 1e-5;
  options->N_beg = 4;
  options->N_end = 16;
  options->N_step = 2;
  options->verbose = 0;
  options->index_function = index_function_all;
}

void
particular_solution_enclosure(geom_t geometry, int angles_coefs[], mpfr_t mu0,
                              mpfr_t nu_low, mpfr_t nu_upp,
                              particular_solution_opt_t options)
{
  mpfr_t *coefs;
  mpfr_t nu, tol, tmp;
  points_t points;
  int prec;

  mpfr_init(nu);
  mpfr_init(tol);
  mpfr_init(tmp);

  prec = mpfr_get_default_prec();

  for (int N = options->N_beg; N <= options->N_end; N += options->N_step)
  {
    /* Compute tolerance and precision to use */
    mpfr_sub(tol, nu_upp, nu_low, MPFR_RNDN);
    mpfr_mul_d(tol, tol, options->tol_relative, MPFR_RNDN);

    /* FIXME: Use the size of nu in the computations for the precision */
    mpfr_log2(tmp, tol, MPFR_RNDN);
    mpfr_mul_d(tmp, tmp, -options->prec_factor, MPFR_RNDN);
    if (mpfr_get_si(tmp, MPFR_RNDN) > prec)
      prec = mpfr_get_si(tmp, MPFR_RNDN);

    mpfr_set_default_prec(prec);

    /* Round variables to new precision */
    mpfr_prec_round(nu_low, prec, MPFR_RNDN);
    mpfr_prec_round(nu_upp, prec, MPFR_RNDN);
    mpfr_prec_round(nu, prec, MPFR_RNDN);

    /* Recompute variables to new precision */
    geom_set(geometry, angles_coefs, prec);

    mpfr_set_prec(mu0, prec);
    mpfr_set_si(mu0, -angles_coefs[1], MPFR_RNDN);
    mpfr_div_si(mu0, mu0, angles_coefs[0], MPFR_RNDN);

    /* Initiate new variables */
    coefs = new mpfr_t[N];
    for (int i = 0; i < N; i++) {
      mpfr_init(coefs[i]);
    }

    points_init(points, 2*N, 2*N);

    boundary(points, geometry, prec);
    interior(points, geometry, prec);

    /* Find the value of nu that minimizes sigma */
    minimize_sigma(nu, points, N, nu_low, nu_upp, mu0, tol,
                   options->index_function);

    mpfr_printf("%i ", N);
    fflush(stdout);

    /* Find the coefficients of the expansion */
    coefs_sigma(coefs, points, N, nu, mu0, options->index_function);

    /* Compute an enclosure of the eigenvalue. To be sure to get correct
       output of the eigenvalue the output of it is handled inside the
       function enclose. */
    enclose(nu_low, nu_upp, angles_coefs, coefs, N, nu,
            options->index_function, 4);
    mpfr_printf("\n");

    for (int i = 0; i < N; i++) {
      mpfr_clear(coefs[i]);
    }
    delete [] coefs;

    points_clear(points);
  }

  mpfr_clear(nu);
  mpfr_clear(tol);
  mpfr_clear(tmp);
}
