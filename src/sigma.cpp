#include "generate-matrix.h"
#include "sigma_eigen.h"

void sigma(mpfr_t res, points_t points, int N, mpfr_t nu, arb_t mu0, int (*index)(int)) {
  mpfr_t *A_mpfr;
  int rows;

  rows = points->boundary + points->interior;

  A_mpfr = new mpfr_t[rows*N];
  for (int i = 0; i < rows*N; i++) {
    mpfr_init(A_mpfr[i]);
  }

  // Fill A_mpfr with the coefficients for the matrix A
  generate_matrix(A_mpfr, points, N, nu, mu0, index);

  sigma_eigen(res, A_mpfr, points->interior, rows, N);

  for (int i = 0; i < rows*N; i++) {
    mpfr_clear(A_mpfr[i]);
  }

  delete [] A_mpfr;
}

void minimize_sigma(mpfr_t nu, points_t points, int N, mpfr_t nu_low,
                    mpfr_t nu_upp, arb_t mu0, mpfr_t tol, int (*index)(int)) {
  mpfr_t a, b, c, d, tmp, tmp2;
  mpfr_t invphi, invphi2, h, yc, yd;
  int n;

  mpfr_init(a);
  mpfr_init(b);
  mpfr_init(c);
  mpfr_init(d);
  mpfr_init(tmp);
  mpfr_init(tmp2);
  mpfr_init(invphi);
  mpfr_init(invphi2);
  mpfr_init(h);
  mpfr_init(yc);
  mpfr_init(yd);

  mpfr_set(a, nu_low, MPFR_RNDN);
  mpfr_set(b, nu_upp, MPFR_RNDN);

  mpfr_sqrt_ui(invphi, 5, MPFR_RNDN);
  mpfr_sub_si(invphi, invphi, 1, MPFR_RNDN);
  mpfr_div_si(invphi, invphi, 2, MPFR_RNDN);

  mpfr_sqrt_ui(invphi2, 5, MPFR_RNDN);
  mpfr_sub_si(invphi2, invphi2, 3, MPFR_RNDN);
  mpfr_div_si(invphi2, invphi2, 2, MPFR_RNDN);
  mpfr_neg(invphi2, invphi2, MPFR_RNDN);

  mpfr_sub(h, b, a, MPFR_RNDN);

  if (mpfr_cmp(h, tol) > 0) {
    mpfr_div(tmp, tol, h, MPFR_RNDN);
    mpfr_log(tmp, tmp, MPFR_RNDN);
    mpfr_log(tmp2, invphi, MPFR_RNDN);
    mpfr_div(tmp, tmp, tmp2, MPFR_RNDN);
    n = mpfr_get_si(tmp, MPFR_RNDU);

    mpfr_mul(tmp, invphi2, h, MPFR_RNDN);
    mpfr_add(c, a, tmp, MPFR_RNDN);
    mpfr_mul(tmp, invphi, h, MPFR_RNDN);
    mpfr_add(d, a, tmp, MPFR_RNDN);

    sigma(yc, points, N, c, mu0, index);
    sigma(yd, points, N, d, mu0, index);

    for (int k = 0; k < n; k++) {
      if (mpfr_cmp(yc, yd) < 0) {
        mpfr_set(b, d, MPFR_RNDN);
        mpfr_set(d, c, MPFR_RNDN);
        mpfr_set(yd, yc, MPFR_RNDN);

        mpfr_mul(h, h, invphi, MPFR_RNDN);

        mpfr_mul(tmp, invphi2, h, MPFR_RNDN);
        mpfr_add(c, a, tmp, MPFR_RNDN);

        sigma(yc, points, N, c, mu0, index);
      } else {
        mpfr_set(a, c, MPFR_RNDN);
        mpfr_set(c, d, MPFR_RNDN);
        mpfr_set(yc, yd, MPFR_RNDN);

        mpfr_mul(h, h, invphi, MPFR_RNDN);

        mpfr_mul(tmp, invphi, h, MPFR_RNDN);
        mpfr_add(d, a, tmp, MPFR_RNDN);

        sigma(yd, points, N, d, mu0, index);
      }
    }

    if (yc < yd) {
      mpfr_set(b, d, MPFR_RNDN);
    } else {
      mpfr_set(a, c, MPFR_RNDN);
      mpfr_set(d, b, MPFR_RNDN);
    }
  }

  mpfr_add(nu, a, b, MPFR_RNDN);
  mpfr_div_si(nu, nu, 2, MPFR_RNDN);

  mpfr_clear(a);
  mpfr_clear(b);
  mpfr_clear(c);
  mpfr_clear(d);
  mpfr_clear(tmp);
  mpfr_clear(tmp2);
  mpfr_clear(invphi);
  mpfr_clear(invphi2);
  mpfr_clear(h);
  mpfr_clear(yc);
  mpfr_clear(yd);
}

void coefs_sigma(mpfr_t *coefs_mpfr, points_t points, int N, mpfr_t nu,
                 arb_t mu0, int (*index)(int)) {
  mpfr_t *A_mpfr;
  int rows;

  rows = points->boundary + points->interior;

  A_mpfr = new mpfr_t[rows*N];
  for (int i = 0; i < rows*N; i++)
  {
    mpfr_init(A_mpfr[i]);
  }

  // Fill A_arr with the coefficients for the matrix A
  generate_matrix(A_mpfr, points, N, nu, mu0, index);

  coefs_sigma_eigen(coefs_mpfr, A_mpfr, points->interior, rows, N);

  for (int i = 0; i < rows*N; i++)
  {
    mpfr_clear(A_mpfr[i]);
  }

  delete [] A_mpfr;
}
