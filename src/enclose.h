#ifndef ENCLOSE
#define ENCLOSE

#include "mpfr.h"

void enclose(arb_t nu_enclosure, int angles_coefs[], arb_ptr coefs,
             int N, arb_t nu, int (*index)(int), int output);

#endif
