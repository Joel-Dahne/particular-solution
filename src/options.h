#ifndef OPTIONS
#define OPTIONS

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
options_struct;

typedef options_struct options_t[1];

void options_default(options_t options);

#endif
