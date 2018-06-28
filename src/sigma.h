#ifndef SIGMA
#define SIGMA

#include "Eigen/Core"
#include "mpreal.h"

using namespace mpfr;

mpreal sigma(mpfr_t *A_arr, struct Points points, mpfr_t *scaling, int N,
             mpfr_t nu, mpfr_t mu0, int (*index)(int));

void coefs_sigma(mpfr_t *coefs_arr, mpfr_t *A_arr, struct Points points,
                 mpfr_t *scaling, int N, mpfr_t nu, mpfr_t mu0, int (*index)(int));

#endif
