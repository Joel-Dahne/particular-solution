#include "generate-matrix.h"
#include "sigma_eigen.h"

void sigma(mpfr_t res, points_t points, slong N, arb_t nu, arb_t mu0,
           int (*index)(int), slong prec) {
  mpfr_t *A_mpfr;
  slong rows;

  rows = points->boundary + points->interior;

  A_mpfr = new mpfr_t[rows*N];
  for (slong i = 0; i < rows*N; i++) {
    mpfr_init(A_mpfr[i]);
  }

  // Fill A_mpfr with the coefficients for the matrix A
  generate_matrix(A_mpfr, points, N, nu, mu0, index, prec);

  sigma_eigen(res, A_mpfr, points->interior, rows, N);

  for (slong i = 0; i < rows*N; i++) {
    mpfr_clear(A_mpfr[i]);
  }

  delete [] A_mpfr;
}

void minimize_sigma(arb_t nu, points_t points, slong N, arb_t nu_enclosure,
                    arb_t mu0, arb_t tol, int (*index)(int), slong prec) {
  mpfr_t yc, yd;
  arb_t a, b, c, d, h, invphi, invphi2, tmp, tmp2;;
  slong n;

  mpfr_init(yc);
  mpfr_init(yd);

  arb_init(a);
  arb_init(b);
  arb_init(c);
  arb_init(d);
  arb_init(invphi);
  arb_init(invphi2);
  arb_init(h);
  arb_init(tmp);
  arb_init(tmp2);

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

    sigma(yc, points, N, c, mu0, index, prec);
    sigma(yd, points, N, d, mu0, index, prec);

    for (slong k = 0; k < n; k++) {
      if (mpfr_cmp(yc, yd) < 0) {
        arb_set(b, d);
        arb_set(d, c);
        mpfr_set(yd, yc, MPFR_RNDN);

        arb_mul(h, h, invphi, prec);

        arb_mul(tmp, invphi2, h, prec);
        arb_add(c, a, tmp, prec);

        sigma(yc, points, N, c, mu0, index, prec);
      } else {
        arb_set(a, c);
        arb_set(c, d);
        mpfr_set(yc, yd, MPFR_RNDN);

        arb_mul(h, h, invphi, prec);

        arb_mul(tmp, invphi, h, prec);
        arb_add(d, a, tmp, prec);

        sigma(yd, points, N, d, mu0, index, prec);
      }
    }

    if (yc < yd) {
      arb_set(b, d);
    } else {
      arb_set(a, c);
      arb_set(d, b);
    }
  }

  arb_add(nu, a, b, prec);
  arb_div_si(nu, nu, 2, prec);

  mpfr_clear(yc);
  mpfr_clear(yd);

  arb_clear(a);
  arb_clear(b);
  arb_clear(c);
  arb_clear(d);
  arb_clear(invphi);
  arb_clear(invphi2);
  arb_clear(h);
  arb_clear(tmp);
  arb_clear(tmp2);
}

void coefs_sigma(arb_ptr coefs, points_t points, slong N, arb_t nu,
                 arb_t mu0, int (*index)(int), slong prec) {
  mpfr_t *A_mpfr, *coefs_mpfr;
  slong rows;

  rows = points->boundary + points->interior;

  A_mpfr = new mpfr_t[rows*N];
  coefs_mpfr = new mpfr_t[N];
  for (slong i = 0; i < rows*N; i++)
  {
    mpfr_init(A_mpfr[i]);
  }
  for (slong i = 0; i < N; i++)
  {
    mpfr_init(coefs_mpfr[i]);
  }

  // Fill A_arr with the coefficients for the matrix A
  generate_matrix(A_mpfr, points, N, nu, mu0, index, prec);

  coefs_sigma_eigen(coefs_mpfr, A_mpfr, points->interior, rows, N);

  for (slong i = 0; i < N; i++)
  {
    arf_set_mpfr(arb_midref(coefs + i), coefs_mpfr[i]);
    mag_zero(arb_radref(coefs + i));
  }

  for (slong i = 0; i < rows*N; i++)
  {
    mpfr_clear(A_mpfr[i]);
  }
  for (slong i = 0; i < N; i++)
  {
    mpfr_clear(coefs_mpfr[i]);
  }

  delete [] A_mpfr;
  delete [] coefs_mpfr;
}
