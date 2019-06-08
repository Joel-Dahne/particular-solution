#ifndef PARTICULAR_SOLUTIONS
#define PARTICULAR_SOLUTIONS

#include "geom.h"
#include "arb.h"

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

void particular_solution_enclosure(arb_t nu_enclosure, geom_t geometry,
                                   particular_solution_opt_t options,
                                   slong prec);

#endif
