#ifndef TOOLS
#define TOOLS

#include "mpfr.h"

void angles_to_vectors(mpfr_t *v1, mpfr_t *v2, mpfr_t *v3, mpfr_t theta_bound,
                       mpfr_t *angles);

void boundary(mpfr_t *thetas, mpfr_t *phis, mpfr_t *v1, mpfr_t *v2, int n,
              int half_boundary);

void interior(mpfr_t *thetas, mpfr_t *phis, mpfr_t *v1, mpfr_t *v2, mpfr_t *v3,
              int n);

#endif
