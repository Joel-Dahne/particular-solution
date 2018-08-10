#include "generate-matrix.h"

#include "unsupported/Eigen/MPRealSupport"
#include "Eigen/QR"
#include "Eigen/SVD"
#include "mpreal.h"

using namespace mpfr;
using namespace Eigen;

typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;

mpreal sigma(points_t points, int N, mpfr_t nu, mpfr_t mu0, int (*index)(int)) {
  mpfr_t *A_mpfr;
  mpreal *A_mpreal;
  MatrixXmp A, Q;
  HouseholderQR<MatrixXmp> qr;
  BDCSVD<MatrixXmp> svd;
  int rows;

  rows = points->boundary + points->interior;

  A_mpfr = new mpfr_t[rows*N];
  A_mpreal = new mpreal[rows*N];

  for (int i = 0; i < rows*N; i++) {
    mpfr_init(A_mpfr[i]);
  }

  // Fill A_mpfr with the coefficients for the matrix A
  generate_matrix(A_mpfr, points, N, nu, mu0, index);

  // Create an array of mpreals with the same values
  for (int i = 0; i < rows*N; i++)
    A_mpreal[i] = mpreal(A_mpfr[i]);

  // Create the matrix A from the coefficients
  A = Map<MatrixXmp>(A_mpreal, rows, N);

  // Perform QR decomposition
  Q = MatrixXmp::Identity(rows, N);
  qr = HouseholderQR<MatrixXmp>(A);
  Q = qr.householderQ() * Q;

  // Find the smallest singular values
  svd = BDCSVD<MatrixXmp>(Q.block(0, 0, points->interior, N));

  for (int i = 0; i < rows*N; i++) {
    mpfr_clear(A_mpfr[i]);
  }

  delete [] A_mpfr;
  delete [] A_mpreal;
  return svd.singularValues()(N-1);
}

void minimize_sigma(mpfr_t nu, points_t points, int N, mpfr_t nu_low,
                    mpfr_t nu_upp, mpfr_t mu0, mpreal tol, int (*index)(int)) {
  mpfr_t a, b, c, d;
  mpreal invphi, invphi2, h, yc, yd;
  int n;

  mpfr_init(a);
  mpfr_init(b);
  mpfr_init(c);
  mpfr_init(d);

  mpfr_set(a, nu_low, MPFR_RNDN);
  mpfr_set(b, nu_upp, MPFR_RNDN);

  invphi = (sqrt(5) - 1)/2;
  invphi2 = (3 - sqrt(5))/2;

  h = mpreal(b) - mpreal(a);

  if (h > tol) {
    n = (int) ceil(log(tol/h)/log(invphi));
    mpfr_add(c, a, (invphi2*h).mpfr_srcptr(), MPFR_RNDN);
    mpfr_add(d, a, (invphi*h).mpfr_srcptr(), MPFR_RNDN);

    yc = sigma(points, N, c, mu0, index);
    yd = sigma(points, N, d, mu0, index);

    for (int k = 0; k < n; k++) {
      if (yc < yd) {
        mpfr_set(b, d, MPFR_RNDN);
        mpfr_set(d, c, MPFR_RNDN);
        yd = yc;
        h = invphi*h;
        mpfr_add(c, a, (invphi2*h).mpfr_srcptr(), MPFR_RNDN);
        yc = sigma(points, N, c, mu0, index);
      } else {
        mpfr_set(a, c, MPFR_RNDN);
        mpfr_set(c, d, MPFR_RNDN);
        yc = yd;
        h = invphi*h;
        mpfr_add(d, a, (invphi*h).mpfr_srcptr(), MPFR_RNDN);
        yd = sigma(points, N, d, mu0, index);
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
}

void coefs_sigma(mpfr_t *coefs_mpfr, points_t points, int N, mpfr_t nu,
                 mpfr_t mu0, int (*index)(int)) {
  mpfr_t *A_mpfr;
  mpreal *A_mpreal;
  MatrixXmp A, Q, coefs;
  HouseholderQR<MatrixXmp> qr;
  BDCSVD<MatrixXmp> svd;
  int rows;

  rows = points->boundary + points->interior;

  A_mpfr = new mpfr_t[rows*N];
  A_mpreal = new mpreal[rows*N];

  for (int i = 0; i < rows*N; i++) {
    mpfr_init(A_mpfr[i]);
  }

  // Fill A_arr with the coefficients for the matrix A
  generate_matrix(A_mpfr, points, N, nu, mu0, index);

  // Create an array of mpreals with the same values
  for (int i = 0; i < rows*N; i++)
    A_mpreal[i] = mpreal(A_mpfr[i]);

  // Create the matrix A from the coefficients
  A = Map<MatrixXmp>(A_mpreal, rows, N);

  // Perform QR decomposition
  Q = MatrixXmp::Identity(rows, N);
  qr = HouseholderQR<MatrixXmp>(A);
  Q = qr.householderQ() * Q;

  // Compute the right singular vector corresponding to the smallest
  // singular value
  svd = BDCSVD<MatrixXmp>(Q.block(0, 0, points->interior, N), ComputeThinV);

  // Compute the coefficients
  coefs = qr.solve(Q*(svd.matrixV().col(N-1)));;

  coefs /= coefs(0);

  for (int i = 0; i < N; i++)
    mpfr_set(coefs_mpfr[i], coefs(i).mpfr_srcptr(), MPFR_RNDN);

  for (int i = 0; i < rows*N; i++) {
    mpfr_clear(A_mpfr[i]);
  }

  delete [] A_mpfr;
  delete [] A_mpreal;
}
