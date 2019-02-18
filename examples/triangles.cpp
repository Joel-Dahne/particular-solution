/*
  Computation of eigenvalues for spherical triangles.

  Author: Joel Dahne
  This file is in the public domain.
*/

#include "particular_solution.h"

#include <string.h>

/* ------------------------------------------------------------------------- */
/*  Different triangles                                                      */
/* ------------------------------------------------------------------------- */

void
get_triangle(geom_t geometry, arb_t nu_enclosure,
             particular_solution_opt_t options, int triangle)
{
  slong angles[6];

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
    arb_set_d(options->prec_factor, 2);

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 0;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 3.056691018);
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
    arb_set_d(options->prec_factor, 2);

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 0;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 3.240902298);
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
    arb_set_d(options->prec_factor, 2);

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 0;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 4.063109028);
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
    arb_set_d(options->prec_factor, 1.2);

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 1;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 4.143210850);
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
    arb_set_d(options->prec_factor, 1.2);

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 0;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 4.470604591);
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
    arb_set_d(options->prec_factor, 1.2);

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 1;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 6.525663100);
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
    arb_set_d(options->prec_factor, 1.2);

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 1;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 1.624084509);
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
    arb_set_d(options->prec_factor, 1.2);

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 1;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 1.825757081);
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
    arb_set_d(options->prec_factor, 1.2);

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 0;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 2.047890892);
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
    arb_set_d(options->prec_factor, 1.2);

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 1;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 2.150869291);
  }

  /* Set up the geometry */
  geom_set_angles(geometry, angles);

  /* Set starting enclosure of nu */
  mag_set_d(arb_radref(nu_enclosure), 0.5);
}

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
  particular_solution_opt_t options;
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

  arb_init(nu_enclosure);

  geom_init(geometry);

  particular_solution_opt_init(options);
  particular_solution_opt_default(options);

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
      arb_set_d(options->tol_relative, atof(argv[i + 1]));
    }
  }

  for (int i = ifrom; i <= ito; i++)
  {
    get_triangle(geometry, nu_enclosure, options, i);

    particular_solution_enclosure(nu_enclosure, geometry, options, prec);
  }

  arb_clear(nu_enclosure);

  geom_clear(geometry);

  particular_solution_opt_clear(options);

  flint_cleanup();

  return 0;
}
