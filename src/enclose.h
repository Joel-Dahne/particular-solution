#ifndef ENCLOSE
#define ENCLOSE

#include "arb.h"
#include "geom.h"

void enclose(arb_t nu_enclosure, geom_t geometry, arb_ptr coefs, int N,
             arb_t nu, int (*index)(int), int output);

#endif
