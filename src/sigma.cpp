#include "generate-matrix.h"

#include "unsupported/Eigen/MPRealSupport"
#include "Eigen/QR"
#include "Eigen/SVD"
#include "mpreal.h"

using namespace mpfr;
using namespace Eigen;

/* Computes the sigma value corresponding to how good of
 an approximation of the eigenvalue that nu is.

 A_arr - Array to store the coefficients of A in, should be allocated
 to have room for (boundary + interior)*N elements.

 thetas, phis - Arrays to contain the theta and phi values to use for
 the points. Should be pointers to Arb vectors. Should have length
 boundary+interior.

 boundary - Number of points corresponding to the boundary.
 interior - Number of points corresponding to the interior.

 N - Number of terms to use in the expansion.

 nu - Value of nu to use in the computations
 */
mpreal sigma(mpfr_t *A_arr, mpfr_t *thetas, mpfr_t *phis, mpfr_t *scaling,
             int boundary, int interior, int N, mpfr_t nu, mpfr_t mu0) {
  typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;
  mpreal A_arr_mpreal[(boundary + interior)*N];

  // Fill A_arr with the coefficients for the matrix A
  generate_matrix(A_arr, thetas, phis, scaling, boundary + interior, N, nu, mu0);

  // Create an array of mpreals with the same values
  for (int i = 0; i < (boundary + interior)*N; i++) {
    A_arr_mpreal[i] = mpreal(A_arr[i]);
  }

  // Create the matrix A from the coefficients
  Map<MatrixXmp> A(A_arr_mpreal, boundary + interior, N);

  // Perform QR decomposition
  MatrixXmp Q(MatrixXmp::Identity(boundary + interior, N));
  HouseholderQR<MatrixXmp> qr(A);
  Q = qr.householderQ() * Q;

  // Find the smallest singular values
  BDCSVD<MatrixXmp> svd(Q.block(0, 0, interior, N));

  return svd.singularValues()(N-1);
}

void coefs_sigma(mpfr_t *coefs_arr, mpfr_t *A_arr, mpfr_t *thetas, mpfr_t *phis,
                 mpfr_t *scaling, int boundary, int interior, int N,
                  mpfr_t nu, mpfr_t mu0) {
  typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;
  mpreal A_arr_mpreal[(boundary + interior)*N];

  // Fill A_arr with the coefficients for the matrix A
  generate_matrix(A_arr, thetas, phis, scaling, boundary + interior, N, nu, mu0);

  // Create an array of mpreals with the same values
  for (int i = 0; i < (boundary + interior)*N; i++) {
    A_arr_mpreal[i] = mpreal(A_arr[i]);
  }

  // Create the matrix A from the coefficients
  Map<MatrixXmp> A(A_arr_mpreal, boundary + interior, N);

  // Perform QR decomposition
  MatrixXmp Q(MatrixXmp::Identity(boundary + interior, N));
  HouseholderQR<MatrixXmp> qr(A);
  Q = qr.householderQ() * Q;

  // Find the smallest singular values
  BDCSVD<MatrixXmp> svd(Q.block(0, 0, interior, N), ComputeThinV);

  MatrixXmp coefs;

  coefs = qr.solve(Q*(svd.matrixV().col(N-1)));

  for (int i = 0; i < N; i++) {
    coefs(i) *= mpreal(scaling[i]);
  }

  coefs /= coefs(0);

  for (int i = 0; i < N; i++) {
    mpfr_set(coefs_arr[i], coefs(i).mpfr_srcptr(), MPFR_RNDN);
  }
}
