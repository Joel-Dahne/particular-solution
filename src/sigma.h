#ifndef SIGMA
#define SIGMA

#include "Eigen/Core"
#include "mpreal.h"

using namespace mpfr;

mpreal sigma(mpfr_t *A_arr, mpfr_t *thetas, mpfr_t *phis, mpfr_t *scaling,
             int boundary, int interior, int N, mpfr_t nu, mpfr_t mu0,
             int (*index)(int));

void coefs_sigma(mpfr_t *coefs_arr, mpfr_t *A_arr, mpfr_t *thetas, mpfr_t *phis,
                 mpfr_t *scaling, int boundary, int interior, int N,
                 mpfr_t nu, mpfr_t mu0, int (*index)(int));

#endif
