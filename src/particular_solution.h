#ifndef PARTICULAR_SOLUTIONS
#define PARTICULAR_SOLUTIONS

#include "geom.h"
#include "arb.h"

typedef struct
{
    arb_t prec_factor;
    arb_t tol_relative;
    slong N_beg;
    slong N_end;
    slong N_step;
    int verbose;
    int (*index_function)(int);

}
particular_solution_opt_struct;

typedef particular_solution_opt_struct particular_solution_opt_t[1];

int index_function_odd(int k);

int index_function_all(int k);

void particular_solution_opt_init(particular_solution_opt_t options);

void particular_solution_opt_default(particular_solution_opt_t options);

void particular_solution_opt_clear(particular_solution_opt_t options);

void particular_solution_enclosure(arb_t nu_enclosure, geom_t geometry,
                                   particular_solution_opt_t options,
                                   slong prec);

#endif
