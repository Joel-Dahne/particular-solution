#ifndef PARTICULAR_SOLUTIONS
#define PARTICULAR_SOLUTIONS

#include "geom.h"
#include "setup.h"
#include "arb.h"

void particular_solution_enclosure(arb_t nu_enclosure, geom_t geometry,
                                   options_t options, slong prec);

#endif
