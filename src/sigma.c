#include "generate_matrix.h"
#include "sigma_eigen.hpp"

void sigma(arb_t res, geom_t geom, points_t points, slong N, arb_t nu,
           slong prec) {
  arb_mat_t A;
  mpfr_t *A_mpfr;
  mpfr_t res_mpfr;
  slong rows, columns;

  mpfr_set_default_prec(prec);

  rows = points->total;
  columns = N*(geom->vertices[0] + geom->vertices[1] + geom->vertices[2]);

  arb_mat_init(A, rows, columns);

  mpfr_init(res_mpfr);

  A_mpfr = (mpfr_t*) calloc(rows*columns, sizeof(mpfr_t));
  for (slong i = 0; i < rows*columns; i++) {
    mpfr_init(A_mpfr[i]);
  }

  /* Fill A with the coefficients for the matrix */
  generate_matrix(A, geom, points, N, nu, prec);

  /* Copy the values to the mpfr version of the matrix */
  for (slong i = 0; i < rows; i++) {
    for (slong j = 0; j < columns; j++) {
      arf_get_mpfr(A_mpfr[j*rows + i], arb_midref(arb_mat_entry(A, i, j)),
                   MPFR_RNDN);
    }
  }

  sigma_eigen(res_mpfr, A_mpfr, points->boundary, rows, columns);

  arf_set_mpfr(arb_midref(res), res_mpfr);

  arb_mat_clear(A);

  mpfr_clear(res_mpfr);

  for (slong i = 0; i < rows*N; i++) {
    mpfr_clear(A_mpfr[i]);
  }

  free(A_mpfr);
}

void minimize_sigma(arb_t nu, geom_t geom, points_t points, slong N,
                    arb_t nu_enclosure, arb_t tol, slong prec) {
  arb_t a, b, c, d, yc, yd, h, invphi, invphi2, tmp, tmp2;;
  slong n;

  arb_init(a);
  arb_init(b);
  arb_init(c);
  arb_init(d);
  arb_init(yc);
  arb_init(yd);
  arb_init(invphi);
  arb_init(invphi2);
  arb_init(h);
  arb_init(tmp);
  arb_init(tmp2);

  n = 0;

  arb_get_lbound_arf(arb_midref(a), nu_enclosure, prec);
  arb_get_ubound_arf(arb_midref(b), nu_enclosure, prec);

  arb_sqrt_ui(invphi, 5, prec);
  arb_sub_si(invphi, invphi, 1, prec);
  arb_div_si(invphi, invphi, 2, prec);

  arb_sqrt_ui(invphi2, 5, prec);
  arb_sub_si(invphi2, invphi2, 3, prec);
  arb_div_si(invphi2, invphi2, 2, prec);
  arb_neg(invphi2, invphi2);

  arb_sub(h, b, a, prec);

  if (arb_ge(h, tol)) {
    arb_div(tmp, tol, h, prec);
    arb_log(tmp, tmp, prec);
    arb_log(tmp2, invphi, prec);
    arb_div(tmp, tmp, tmp2, prec);
    n = arf_get_si(arb_midref(tmp),  ARF_RND_CEIL);

    arb_mul(tmp, invphi2, h, prec);
    arb_add(c, a, tmp, prec);
    arb_mul(tmp, invphi, h, prec);
    arb_add(d, a, tmp, prec);

    sigma(yc, geom, points, N, c, prec);
    sigma(yd, geom, points, N, d, prec);

    for (slong k = 0; k < n; k++) {
      if (arf_cmp(arb_midref(yc), arb_midref(yd)) < 0) {
        arb_set(b, d);
        arb_set(d, c);
        arb_set(yd, yc);

        arb_mul(h, h, invphi, prec);

        arb_mul(tmp, invphi2, h, prec);
        arb_add(c, a, tmp, prec);

        sigma(yc, geom, points, N, c, prec);
      } else {
        arb_set(a, c);
        arb_set(c, d);
        arb_set(yc, yd);

        arb_mul(h, h, invphi, prec);

        arb_mul(tmp, invphi, h, prec);
        arb_add(d, a, tmp, prec);

        sigma(yd, geom, points, N, d, prec);
      }
    }

    if (arf_cmp(arb_midref(yc), arb_midref(yd)) < 0) {
      arb_set(b, d);
    } else {
      arb_set(a, c);
    }
  }

  arb_add(nu, a, b, prec);
  arb_div_si(nu, nu, 2, prec);

#ifdef DEBUG
  flint_printf("DEBUG minimize_sigma: Iterations %d, minimum: ", n);
  arf_printd(arb_midref(nu), (slong)ceil(prec*log10(2)));
  flint_printf("\n");
#endif

  arb_clear(a);
  arb_clear(b);
  arb_clear(c);
  arb_clear(d);
  arb_clear(yc);
  arb_clear(yd);
  arb_clear(invphi);
  arb_clear(invphi2);
  arb_clear(h);
  arb_clear(tmp);
  arb_clear(tmp2);
}

void coefs_sigma(arb_ptr* coefs, geom_t geom, points_t points, slong N,
                 arb_t nu, slong prec) {
  arb_mat_t A;
  mpfr_t *A_mpfr, *coefs_mpfr;
  slong rows, columns;
  slong start;

  mpfr_set_default_prec(prec);

  rows = points->total;
  columns = N*(geom->vertices[0] + geom->vertices[1] + geom->vertices[2]);

  arb_mat_init(A, rows, columns);

  A_mpfr = (mpfr_t*) calloc(rows*columns, sizeof(mpfr_t));
  coefs_mpfr = (mpfr_t*) calloc(columns, sizeof(mpfr_t));
  for (slong i = 0; i < rows*columns; i++)
  {
    mpfr_init(A_mpfr[i]);
  }
  for (slong i = 0; i < columns; i++)
  {
    mpfr_init(coefs_mpfr[i]);
  }

  /* Fill A with the coefficients for the matrix */
  generate_matrix(A, geom, points, N, nu, prec);

  /* Copy the values to the mpfr version of the matrix */
  for (slong i = 0; i < rows; i++) {
    for (slong j = 0; j < columns; j++) {
      arf_get_mpfr(A_mpfr[j*rows + i], arb_midref(arb_mat_entry(A, i, j)),
                   MPFR_RNDN);
    }
  }

  coefs_sigma_eigen(coefs_mpfr, A_mpfr, points->boundary, rows, columns);

  start = 0;
  for (slong vertex = 0; vertex < 3; vertex++)
  {
    if (geom->vertices[vertex])
    {
      for (slong i = 0; i < N; i++)
      {
        arf_set_mpfr(arb_midref(coefs[vertex] + i), coefs_mpfr[start + i]);
        mag_zero(arb_radref(coefs[vertex] + i));
      }

      start += N;
    }
  }

  arb_mat_clear(A);

  for (slong i = 0; i < rows*columns; i++)
  {
    mpfr_clear(A_mpfr[i]);
  }
  for (slong i = 0; i < columns; i++)
  {
    mpfr_clear(coefs_mpfr[i]);
  }

  free(A_mpfr);
  free(coefs_mpfr);
}
