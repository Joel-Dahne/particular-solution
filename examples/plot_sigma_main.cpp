/*
  Plotting sigma values on an interval computed when using the method
  of particular solutions on spherical triangles.

  Author: Joel Dahne
  This file is in the public domain.
*/

#include "options.h"
#include "setup.h"
#include "plot_sigma.h"

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
  arb_t tmp;
  arf_t inf, sup;
  geom_t geometry;
  options_t options;
  slong triangle, num_points, prec;

  srand(1);

  triangle = -1;

  for (slong i = 1; i < argc; i++)
  {
    if (!strcmp(argv[i], "-i"))
    {
      triangle = atol(argv[i+1]);
      if (triangle < 0 || triangle >= NUM_TRIANGLES)
        flint_abort();
    }
  }

  if (triangle == -1)
  {
    mpfr_printf("Plot sigma values on an interval using plot_sigma.\n");
    mpfr_printf("Usage: plot_sigma_main -i n [-prec p] [-N_beg b] [...]\n\n");
    mpfr_printf("-i n       - compute for triangle n (0 <= n <= %d)\n", NUM_TRIANGLES - 1);
    mpfr_printf("-inf a    - lower bound for interval (default 0.0)\n");
    mpfr_printf("-sup b    - upper bound for interval (default 10.0)\n");
    mpfr_printf("-prec p    - precision in bits (default p = 64)\n");
    mpfr_printf("-N_beg b    - N value to start at (default 4)\n");
    mpfr_printf("-N_end e    - N value to stop at (default 16)\n");
    mpfr_printf("-N_step s    - step to take with N each iteration (default 2)\n");
    mpfr_printf("-points num  - number of points to use in the plot (default 500)\n");
    mpfr_printf("Implemented triangles:\n");
    for (int i = 0; i < NUM_TRIANGLES; i++)
      mpfr_printf("T%d = %s\n", i, descr[i]);
    mpfr_printf("\n");
    return 1;
  }

  num_points = 500;
  prec = 64;

  arb_init(tmp);

  arf_init(inf);
  arf_init(sup);

  geom_init(geometry);

  options_init(options);
  options_default(options);

  arf_set_si(inf, 0);
  arf_set_si(sup, 10);

  for (int i = 1; i < argc; i++)
  {
    if (!strcmp(argv[i], "-inf"))
    {
      arf_set_d(inf, atof(argv[i + 1]));
    }
    else if (!strcmp(argv[i], "-sup"))
    {
      arf_set_d(sup, atof(argv[i + 1]));
    }
    else if (!strcmp(argv[i], "-prec"))
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
    else if (!strcmp(argv[i], "-points"))
    {
      num_points = atol(argv[i + 1]);
    }
  }

  get_domain(geometry, tmp, options, triangle);
  plot_sigma(inf, sup, geometry, num_points, options, prec);

  arb_clear(tmp);

  arf_clear(inf);
  arf_clear(sup);

  geom_clear(geometry);

  options_clear(options);

  flint_cleanup();

  return 0;
}
