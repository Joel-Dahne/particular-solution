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
  options->verbose = 0;
  options->index_function = index_function_all;
}

void
particular_solution_opt_clear(particular_solution_opt_t options)
{
  arb_clear(options->prec_factor);
  arb_clear(options->tol_relative);
}

void
particular_solution_enclosure(geom_t geometry, arb_t nu_enclosure,
                              particular_solution_opt_t options)
{
  arb_ptr coefs;
  arb_t nu, mu0, tol, tmp;
  points_t points;
  int prec;

  arb_init(nu);
  arb_init(mu0);
  arb_init(tol);
  arb_init(tmp);

  prec = mpfr_get_default_prec();

  for (int N = options->N_beg; N <= options->N_end; N += options->N_step)
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

    arb_set_fmpq(mu0, geometry->angles, prec);
    arb_inv(mu0, mu0, prec);
    arb_neg(mu0, mu0);

    /* Initiate new variables */
    coefs = _arb_vec_init(N);

    points_init(points, 2*N, 2*N);

    boundary(points, geometry, prec);
    interior(points, geometry, prec);

    /* Find the value of nu that minimizes sigma */
    minimize_sigma(nu, points, N, nu_enclosure, mu0, tol,
                   options->index_function);

    flint_printf("%i ", N);
    fflush(stdout);

    /* Find the coefficients of the expansion */
    coefs_sigma(coefs, points, N, nu, mu0, options->index_function);

    /* Compute an enclosure of the eigenvalue. To be sure to get correct
       output of the eigenvalue the output of it is handled inside the
       function enclose. */
    enclose(nu_enclosure, geometry, coefs, N, nu, options->index_function,
            4);
    flint_printf("\n");

    _arb_vec_clear(coefs, N);

    points_clear(points);
  }

  arb_clear(nu);
  arb_clear(mu0);
  arb_clear(tol);
  arb_clear(tmp);
}
