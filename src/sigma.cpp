#include "generate-matrix.h"

#include "unsupported/Eigen/MPRealSupport"
#include "Eigen/QR"
#include "Eigen/SVD"
#include "mpreal.h"

using namespace mpfr;
using namespace Eigen;

mpreal sigma(mpfr_t *A_arr, points_t points, int N, mpfr_t nu, mpfr_t mu0,
             int (*index)(int)) {
  typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;
  mpreal *A_arr_mpreal;
  A_arr_mpreal = new mpreal[(points->boundary + points->interior)*N];

  // Fill A_arr with the coefficients for the matrix A
  generate_matrix(A_arr, points, N, nu, mu0, index);

  // Create an array of mpreals with the same values
  for (int i = 0; i < (points->boundary + points->interior)*N; i++)
    A_arr_mpreal[i] = mpreal(A_arr[i]);

  // Create the matrix A from the coefficients
  Map<MatrixXmp> A(A_arr_mpreal, points->boundary + points->interior, N);

  // Perform QR decomposition
  MatrixXmp Q(MatrixXmp::Identity(points->boundary + points->interior, N));
  HouseholderQR<MatrixXmp> qr(A);
  Q = qr.householderQ() * Q;

  // Find the smallest singular values
  BDCSVD<MatrixXmp> svd(Q.block(0, 0, points->interior, N));

  delete [] A_arr_mpreal;
  return svd.singularValues()(N-1);
}

void minimize_sigma(mpfr_t nu, mpfr_t *A_arr, points_t points, int N,
                    mpfr_t nu_low, mpfr_t nu_upp, mpfr_t mu0, mpreal tol,
                    int (*index)(int)) {
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

    yc = sigma(A_arr, points, N, c, mu0, index);
    yd = sigma(A_arr, points, N, d, mu0, index);

    for (int k = 0; k < n; k++) {
      if (yc < yd) {
        mpfr_set(b, d, MPFR_RNDN);
        mpfr_set(d, c, MPFR_RNDN);
        yd = yc;
        h = invphi*h;
        mpfr_add(c, a, (invphi2*h).mpfr_srcptr(), MPFR_RNDN);
        yc = sigma(A_arr, points, N, c, mu0, index);
      } else {
        mpfr_set(a, c, MPFR_RNDN);
        mpfr_set(c, d, MPFR_RNDN);
        yc = yd;
        h = invphi*h;
        mpfr_add(d, a, (invphi*h).mpfr_srcptr(), MPFR_RNDN);
        yd = sigma(A_arr, points, N, d, mu0, index);
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

void coefs_sigma(mpfr_t *coefs_arr, mpfr_t *A_arr, points_t points,
                 int N, mpfr_t nu, mpfr_t mu0, int (*index)(int)) {
  typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;
  mpreal *A_arr_mpreal;
  A_arr_mpreal = new mpreal[(points->boundary + points->interior)*N];

  // Fill A_arr with the coefficients for the matrix A
  generate_matrix(A_arr, points, N, nu, mu0, index);

  // Create an array of mpreals with the same values
  for (int i = 0; i < (points->boundary + points->interior)*N; i++)
    A_arr_mpreal[i] = mpreal(A_arr[i]);

  // Create the matrix A from the coefficients
  Map<MatrixXmp> A(A_arr_mpreal, points->boundary + points->interior, N);

  // Perform QR decomposition
  MatrixXmp Q(MatrixXmp::Identity(points->boundary + points->interior, N));
  HouseholderQR<MatrixXmp> qr(A);
  Q = qr.householderQ() * Q;

  // Compute the right singular vector corresponding to the smallest
  // singular value
  BDCSVD<MatrixXmp> svd(Q.block(0, 0, points->interior, N), ComputeThinV);

  // Compute the coefficients
  MatrixXmp coefs = qr.solve(Q*(svd.matrixV().col(N-1)));;

  coefs /= coefs(0);

  for (int i = 0; i < N; i++)
    mpfr_set(coefs_arr[i], coefs(i).mpfr_srcptr(), MPFR_RNDN);

  delete [] A_arr_mpreal;
}
