#ifndef ENCLOSE
#define ENCLOSE

#include "mpfr.h"

void enclose(mpfr_t nu_low, mpfr_t nu_upp, int angles_coefs[], arb_ptr coefs,
             int N, arb_t nu, int (*index)(int), int output);

#endif
