#include "enclose.h"

#include "norm.h"
#include "eigenfunction.h"
#include "maximum.h"

void
enclose(arb_t nu_enclosure, geom_t geom, arb_ptr* coefs,
        slong N, arb_t nu, slong vertex, slong prec) {
  arb_t eps, norm, max, eigenvalue, tmp;

  arb_init(eps);
  arb_init(norm);
  arb_init(max);
  arb_init(eigenvalue);
  arb_init(tmp);

  /* Compute an enclosure of a lower bound of the norm */
  integral_norm(norm, geom, coefs[vertex], N, nu, vertex, prec);

  /* Compute an enclosure of the maximum on the boundary */
  maximize(max, geom, coefs[vertex], N, nu, vertex, prec);

  /* Set eps equal to the square root of the area of the triangle */
  geom_area(eps, geom, prec);
  arb_sqrt(eps, eps, prec);

  /* Multiply with the maximum and divide with the norm */
  arb_mul(eps, eps, max, prec);
  arb_div(eps, eps, norm, prec);

  /* Compute lower and upper bounds for the eigenvalue */
  arb_add_si(eigenvalue, nu, 1, prec);
  arb_mul(eigenvalue, eigenvalue, nu, prec);
  arb_set_si(tmp, 1);
  arb_add_error(tmp, eps);
  arb_div(eigenvalue, eigenvalue, tmp, prec);

  /* Compute lower and upper bounds for the enclosure of nu */
  arb_set_d(tmp, 0.25);
  arb_add(tmp, tmp, eigenvalue, prec);
  arb_sqrt(tmp, tmp, prec);
  arb_set_d(eps, 0.5);
  arb_sub(tmp, tmp, eps, prec);

  if (arb_is_finite(tmp) && arb_overlaps(nu_enclosure, tmp))
  {
    arb_intersection(nu_enclosure, nu_enclosure, tmp, prec);
  }

  arb_clear(eps);
  arb_clear(norm);
  arb_clear(max);
  arb_clear(eigenvalue);
  arb_clear(tmp);
}
