#include "options.h"

void
options_default(options_t options)
{
  options->tol_relative = 1e-5;
  options->prec_factor = 1.2;
  options->N_beg = 4;
  options->N_end = 16;
  options->N_step = 2;
  options->plot_n = 100;
  options->output = 0;
  options->output_final = 0;
  options->output_time = 0;
}
