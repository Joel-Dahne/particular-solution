/*
  Computation of eigenvalues for spherical triangles.

  Author: Joel Dahne
  This file is in the public domain.
*/

#include "options.h"
#include "setup.h"
#include "particular_solution.h"

#include <string.h>

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

/* ------------------------------------------------------------------------- */
/*  Main test program                                                        */
/* ------------------------------------------------------------------------- */

int
main(int argc, char *argv[])
{
  arb_t nu_enclosure;
  geom_t geometry;
  options_t options;
  slong ifrom, ito, prec;

  ifrom = ito = -1;

  for (slong i = 1; i < argc; i++)
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
    mpfr_printf("-output m    - output type\n");
    mpfr_printf("             - 0: no output\n");
    mpfr_printf("             - 1: enclosure for eigenvalue\n");
    mpfr_printf("             - 2: enclosure for nu\n");
    mpfr_printf("             - 3: midpoint of enclosure for eigenvalue\n");
    mpfr_printf("             - 4: midpoint of enclosure for nu\n");
    mpfr_printf("             - 5: width of enclosure for eigenvalue\n");
    mpfr_printf("             - 6: width of enclosure for nu\n");
    mpfr_printf("             - 7: coefficients for the approximate eigenfunction\n");
    mpfr_printf("             - 8: plot of the approximate eigenfunction\n");
    mpfr_printf("             - 9: minimizing nu value\n");
    mpfr_printf("-prec p      - precision in bits (default p = 64)\n");
    mpfr_printf("-N_beg b     - N value to start at (default 4)\n");
    mpfr_printf("-N_end e     - N value to stop at (default 16)\n");
    mpfr_printf("-N_step s    - step to take with N each iteration (default 2)\n");
    mpfr_printf("-tol eps     - relative accuracy goal each iteration (default 1e-5)\n");
    mpfr_printf("-plot n      - determines how to enclose eigenfunction with output 8 (default 0)\n");
    mpfr_printf("             -  -1: non-rigorous plot on points on the boundary\n");
    mpfr_printf("             -   0: simple interval enclosure\n");
    mpfr_printf("             - > 0: enclosure with Taylor expansion with n terms\n");
    mpfr_printf("-final       - flag for outputting only for the last N value (default off)\n");
    mpfr_printf("-time       - flag for outputting timing information (default off)\n");
    mpfr_printf("Implemented triangles:\n");
    for (int i = 0; i < NUM_TRIANGLES; i++)
      mpfr_printf("T%d = %s\n", i, descr[i]);
    mpfr_printf("\n");
    return 1;
  }

  prec = 64;

  arb_init(nu_enclosure);

  geom_init(geometry);

  options_default(options);

  options->output = 1;

  for (int i = 1; i < argc; i++)
  {
    if (!strcmp(argv[i], "-prec"))
    {
      prec = atol(argv[i + 1]);
    }
    else if (!strcmp(argv[i], "-output"))
    {
      options->output = atol(argv[i + 1]);
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
    else if (!strcmp(argv[i], "-plot"))
    {
      options->plot_n = atol(argv[i + 1]);
    }
    else if (!strcmp(argv[i], "-final"))
    {
      options->output_final = 1;
    }
    else if (!strcmp(argv[i], "-time"))
    {
      options->output_time = 1;
    }
  }

  for (int i = ifrom; i <= ito; i++)
  {
    get_domain(geometry, nu_enclosure, options, i);

    particular_solution_enclosure(nu_enclosure, geometry, options, prec);
  }

  arb_clear(nu_enclosure);

  geom_clear(geometry);

  flint_cleanup();

  return 0;
}
