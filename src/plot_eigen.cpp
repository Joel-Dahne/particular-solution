#include "plot_eigen.h"

#include "generate-matrix.h"

void
plot_eigen(geom_t geometry, arb_ptr coefs, slong N, arb_t nu, slong num_points,
           slong prec) {
  arb_ptr evals;
  points_t  points;

  evals = _arb_vec_init(num_points);

  points_init(points, num_points, 0);
  boundary(points, geometry, prec);

  eigenfunction(evals, geometry, coefs, points, N, nu, prec);

  for (slong i = 0; i < num_points; i++)
  {
    arf_printd(arb_midref(evals + i), (slong)ceil(prec*log10(2)));
    flint_printf("\n");
  }

  _arb_vec_clear(evals, num_points);

  points_clear(points);
}
