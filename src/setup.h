#ifndef SETUP
#define SETUP

#include "geom.h"
#include "arb.h"

/* Stores options used in the computations. */
typedef struct
{
  double tol_relative;
  double prec_factor;
  slong N_beg;
  slong N_end;
  slong N_step;
  slong plot_n;
  int output;
  int output_final;
  int output_time;
}
particular_solution_opt_struct;

typedef particular_solution_opt_struct particular_solution_opt_t[1];

void particular_solution_opt_default(particular_solution_opt_t options);

void get_triangle_defaults(geom_t geometry, arb_t nu_enclosure,
                           particular_solution_opt_t options, int triangle);

#endif
