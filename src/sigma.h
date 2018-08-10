#ifndef SIGMA
#define SIGMA

#include "Eigen/Core"
#include "mpreal.h"

using namespace mpfr;

mpreal sigma(points_t points, int N, mpfr_t nu, mpfr_t mu0, int (*index)(int));

void minimize_sigma(mpfr_t nu, points_t points, int N, mpfr_t nu_low,
                    mpfr_t nu_upp, mpfr_t mu0, mpreal tol, int (*index)(int));

void coefs_sigma(mpfr_t *coefs_mpfr, points_t points, int N, mpfr_t nu,
                 mpfr_t mu0, int (*index)(int));

#endif
