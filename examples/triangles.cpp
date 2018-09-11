/*
  Computation of eigenvalues for spherical triangles with at most one
  singular corner.

  Author: Joel Dahne
  This file is in the public domain.
*/

#include "particular_solution.h"
#include "arb.h"

#include <string.h>

/* ------------------------------------------------------------------------- */
/*  Different triangles                                                      */
/* ------------------------------------------------------------------------- */

void
get_triangle(geom_t geometry, int angles[], mpfr_t mu0, mpfr_t nu_low,
             mpfr_t nu_upp, particular_solution_opt_t options, int triangle,
             int prec)
{
  /* Set default precision */
  mpfr_set_default_prec(prec);

  if (triangle == 0)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 3;
    angles[1] = 4;
    angles[2] = 1;
    angles[3] = 3;
    angles[4] = 1;
    angles[5] = 2;

    /* Set custom options */
    options->prec_factor = 2;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 0;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 3.056691018, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 1)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 2;
    angles[1] = 3;
    angles[2] = 1;
    angles[3] = 3;
    angles[4] = 1;
    angles[5] = 2;

    /* Set custom options */
    options->prec_factor = 2;
    options->index_function = index_function_all;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 0;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 3.240902298, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 2)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 2;
    angles[1] = 3;
    angles[2] = 1;
    angles[3] = 4;
    angles[4] = 1;
    angles[5] = 2;

    /* Set custom options */
    options->prec_factor = 2;
    options->index_function = index_function_all;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 0;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 4.063109028, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 3)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 2;
    angles[1] = 3;
    angles[2] = 1;
    angles[3] = 3;
    angles[4] = 1;
    angles[5] = 3;

    /* Set custom options */
    options->prec_factor = 1.2;
    options->index_function = index_function_all;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 1;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 4.143210850, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 4)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 3;
    angles[1] = 4;
    angles[2] = 1;
    angles[3] = 4;
    angles[4] = 1;
    angles[5] = 3;

    /* Set custom options */
    options->prec_factor = 1.2;
    options->index_function = index_function_all;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 0;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 4.470604591, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 5)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 2;
    angles[1] = 3;
    angles[2] = 1;
    angles[3] = 4;
    angles[4] = 1;
    angles[5] = 4;

    /* Set custom options */
    options->prec_factor = 1.2;
    options->index_function = index_function_odd;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 1;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 6.525663100, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 6)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 2;
    angles[1] = 3;
    angles[2] = 3;
    angles[3] = 4;
    angles[4] = 3;
    angles[5] = 4;

    /* Set custom options */
    options->prec_factor = 1.2;
    options->index_function = index_function_odd;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 1;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 1.624084509, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 7)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 2;
    angles[1] = 3;
    angles[2] = 2;
    angles[3] = 3;
    angles[4] = 2;
    angles[5] = 3;

    /* Set custom options */
    options->prec_factor = 1.2;
    options->index_function = index_function_odd;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 1;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 1.825757081, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 8)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 1;
    angles[1] = 2;
    angles[2] = 2;
    angles[3] = 3;
    angles[4] = 3;
    angles[5] = 4;

    /* Set custom options */
    options->prec_factor = 1.2;
    options->index_function = index_function_all;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 0;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 2.047890892, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 9)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 1;
    angles[1] = 2;
    angles[2] = 2;
    angles[3] = 3;
    angles[4] = 2;
    angles[5] = 3;

    /* Set custom options */
    options->prec_factor = 1.2;
    options->index_function = index_function_odd;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 1;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 2.150869291, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }

  /* Set up the geometry */
  geom_set_prec(geometry, prec);
  geom_set(geometry, angles);

  /* Set mu0 */
  mpfr_set_prec(mu0, prec);
  mpfr_set_si(mu0, -angles[1], MPFR_RNDN);
  mpfr_div_si(mu0, mu0, angles[0], MPFR_RNDN);

  /* Set starting enclosure of nu */
  mpfr_sub_d(nu_low, nu_low, 1e-2, MPFR_RNDN);
  mpfr_add_d(nu_upp, nu_upp, 1e-2, MPFR_RNDN);
}

/* ------------------------------------------------------------------------- */
/*  Main test program                                                        */
/* ------------------------------------------------------------------------- */

#define NUM_TRIANGLES 10

const char * descr[NUM_TRIANGLES] =
{
  "(3pi/4, pi/3, pi/2)",
  "(2pi/3, pi/3, pi/2)",
  "(2pi/3, pi/4, pi/2)",
  "(2pi/3, pi/3, pi/3)",
  "(3pi/4, pi/4, pi/3)",
  "(2pi/3, pi/4, pi/4)",
  "(2pi/3, 3pi/4, 3pi/4)",
  "(2pi/3, 2pi/3, 2pi/3)",
  "(pi/2, 2pi/3, 3pi/4)",
  "(pi/2, 2pi/3, 2pi/3)"
};

int
main(int argc, char *argv[])
{
  mpfr_t mu0, nu_low, nu_upp;
  int angles[6];
  geom_t geometry;
  particular_solution_opt_t options;
  int ifrom, ito, prec;

  ifrom = ito = -1;

  for (int i = 1; i < argc; i++)
  {
    if (!strcmp(argv[i], "-i"))
    {
      if (!strcmp(argv[i+1], "all"))
      {
        ifrom = 0;
        ito = NUM_TRIANGLES - 1;
      }
      else if (!strcmp(argv[i+1], "regular"))
      {
        ifrom = 0;
        ito = 5;
      }
      else if (!strcmp(argv[i+1], "singular"))
      {
        ifrom = 6;
        ito = 9;
      }
      else
      {
        ifrom = ito = atol(argv[i+1]);
        if (ito < 0 || ito >= NUM_TRIANGLES)
          flint_abort();
      }
    }
  }

  if (ifrom == -1)
  {
    mpfr_printf("Compute eigenvalues using particular_solution.\n");
    mpfr_printf("Usage: triangles -i n [-prec p] [-N_beg b] [...]\n\n");
    mpfr_printf("-i n       - compute for triangle n (0 <= n <= %d), or \"-i all\",\n", NUM_TRIANGLES - 1);
    mpfr_printf("           - or \"-i regular\" or \"-i singular\"\n");
    mpfr_printf("-prec p    - precision in bits (default p = 64)\n");
    mpfr_printf("-N_beg b    - N value to start at (default 4)\n");
    mpfr_printf("-N_end e    - N value to stop at (default 16)\n");
    mpfr_printf("-N_step s    - step to take with N each iteration (default 2)\n");
    mpfr_printf("-tol eps    - relative accuracy goal each iteration (default 1e-5)\n");
    mpfr_printf("Implemented triangles:\n");
    for (int i = 0; i < NUM_TRIANGLES; i++)
      mpfr_printf("T%d = %s\n", i, descr[i]);
    mpfr_printf("\n");
    return 1;
  }

  prec = 64;

  mpfr_init(mu0);
  mpfr_init(nu_low);
  mpfr_init(nu_upp);

  geom_init(geometry);

  particular_solution_opt_init(options);

  for (int i = 1; i < argc; i++)
  {
    if (!strcmp(argv[i], "-prec"))
    {
      prec = atol(argv[i + 1]);
    }
    else if (!strcmp(argv[i], "-N_beg"))
    {
      options->N_beg = atol(argv[i + 1]);
    }
    else if (!strcmp(argv[i], "-N_end"))
    {
      options->N_end = atol(argv[i + 1]);
    }
    else if (!strcmp(argv[i], "-N_step"))
    {
      options->N_step = atol(argv[i + 1]);
    }
    else if (!strcmp(argv[i], "-tol"))
    {
      options->tol_relative = atof(argv[i + 1]);
    }
  }

  for (int i = ifrom; i <= ito; i++)
  {
    get_triangle(geometry, angles, mu0, nu_low, nu_upp, options, i, prec);

    particular_solution_enclosure(geometry, angles, mu0, nu_low, nu_upp,
                                  options);
  }

  mpfr_clear(mu0);
  mpfr_clear(nu_low);
  mpfr_clear(nu_upp);

  geom_clear(geometry);

  flint_cleanup();

  return 0;
}
