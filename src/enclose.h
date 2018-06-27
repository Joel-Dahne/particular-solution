#ifndef ENCLOSE
#define ENCLOSE

#include "mpfr.h"

void enclose(mpfr_t eps, int angles_coefs[], mpfr_t *coefs, int N,
             mpfr_t nu, int index_step);

#endif
