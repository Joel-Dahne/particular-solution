#ifndef PLOT_EIGEN
#define PLOT_EIGEN

#include "geom.h"

void plot_eigen(geom_t geometry, arb_ptr coefs, slong N, arb_t nu, arb_t mu0,
                slong num_points, int (*index)(int), slong prec);


#endif
