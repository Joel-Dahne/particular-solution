#ifndef PARTICULAR_SOLUTIONS
#define PARTICULAR_SOLUTIONS

#include "geom.h"
#include "arb.h"

typedef struct
{
    double prec_factor;
    double tol_relative;
    int N_beg;
    int N_end;
    int N_step;
    int verbose;
    int (*index_function)(int);

}
particular_solution_opt_struct;

typedef particular_solution_opt_struct particular_solution_opt_t[1];

int index_function_odd(int k);

int index_function_all(int k);

void particular_solution_opt_init(particular_solution_opt_t options);

void particular_solution_enclosure(geom_t geometry, int angles_coefs[],
                                   mpfr_t nu_low, mpfr_t nu_upp,
                                   particular_solution_opt_t options);

#endif
