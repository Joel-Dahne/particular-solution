#ifndef PLOT_SIGMA
#define PLOT_SIGMA

#include "geom.h"
#include "particular_solution.h"

void plot_sigma(arf_t inf, arf_t sup, geom_t geometry, slong num_points,
                particular_solution_opt_t options, slong prec);

#endif
