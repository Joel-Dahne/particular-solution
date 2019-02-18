#include "generate-matrix.h"

#include "arb_hypgeom.h"

static void
eval_eigenfunction(arb_t res, arb_t theta, arb_t phi, arb_t nu, arb_t mu,
                   slong prec)
{
  arb_t tmp;
  slong prec_local;

  arb_init(tmp);

  arb_cos(tmp, theta, prec);
  prec_local = prec;
  do {
    arb_hypgeom_legendre_p(res, nu, mu, tmp, 0, prec_local);
    prec_local *=2;
  } while (!arb_is_finite(res));

  arb_mul(tmp, mu, phi, prec);
  arb_sin(tmp, tmp, prec);

  arb_mul(res, res, tmp, prec);

  arb_clear(tmp);
}

/* Generates a matrix for computations of the sigma values. The rows
 * of the matrix corresponds to points in the points argument and the
 * columns corresponds to different functions in the expansion. The
 * value of an element is given by the corresponding function
 * evaluated at the corresponding point.
 *
 * Expansions from up to all three vertices can be used. For each
 * vertex N terms are used in the expansion. The first N columns
 * correspond to the first used vertex, the second N columns to the
 * second used vertex and similar for the third used vertex. The
 * number of columns is therefore equal to N times the number of used
 * vertices. */
void generate_matrix(mpfr_t *A, geom_t geom, points_t points, slong N, arb_t nu,
                     slong prec)
{
  arb_t mu, res;
  slong n, column_start;

  arb_init(mu);
  arb_init(res);

  n = points->total;
  column_start = 0;

  for (slong vertex = 0; vertex < 3; vertex++)
    {
    /* Check if an expansion from the vertex is to be used */
    if (geom->vertices[vertex]) {

      for (slong i = 0; i < n; i++) {
        /* If the values for points->thetas + i and points->phis + i
           are indeterminate set function value to zero. This is what
           points_t expects. */
        if (arb_is_finite(points->thetas[vertex] + i)
            && arb_is_finite(points->phis[vertex] + i)) {

          for (slong j = column_start; j < column_start + N; j++) {
            geom_get_mu(mu, geom, vertex, j - column_start, prec);
            eval_eigenfunction(res, points->thetas[vertex] + i,
                               points->phis[vertex] + i, nu, mu, prec);

            arf_get_mpfr(A[j*n + i], arb_midref(res), MPFR_RNDN);
          }
        } else {
          for (slong j = column_start; j < column_start + N; j++) {
            mpfr_set_zero(A[j*n + i], 0);
          }
        }
      }

      column_start += N;
    }
  }

  arb_clear(mu);
  arb_clear(res);
}

/*   Evaluates the eigenfunction given by the sum of the basis
 *   eigenfunctions with the supplied coefficients at the supplied
 *   points. */
void eigenfunction(arb_ptr res, geom_t geom, arb_ptr coefs, points_t points,
                   slong N, arb_t nu, slong prec) {
  arb_t mu, term;
  slong n;

  arb_init(mu);
  arb_init(term);

  n = points->total;

  for (slong i = 0; i < n; i++) {
    arb_zero(res + i);

    for (slong vertex = 0; vertex < 3; vertex++) {
      if (geom->vertices[vertex]) {
        for (slong j = 0; j < N; j++) {
          /* If the values for points->thetas + i and points->phis + i
             are indeterminate set function value to zero. This is what
             points_t expects. */
          if (arb_is_finite(points->thetas[vertex] + i)
              && arb_is_finite(points->phis[vertex] + i)) {
            geom_get_mu(mu, geom, vertex, j, prec);
            eval_eigenfunction(term, points->thetas[vertex] + i,
                               points->phis[vertex] + i, nu, mu, prec);
            arb_addmul(res + i, term, coefs + j, prec);
          }
        }
      }
    }
  }

  arb_clear(mu);
  arb_clear(term);
}
