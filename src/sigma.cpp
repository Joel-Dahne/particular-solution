#include "generate-matrix.h"

#include "unsupported/Eigen/MPRealSupport"
#include "Eigen/QR"
#include "Eigen/SVD"
#include "mpreal.h"

using namespace mpfr;
using namespace Eigen;

struct Points {
  mpfr_t *thetas;
  mpfr_t *phis;
  int boundary;
  int interior;
};

mpreal sigma(mpfr_t *A_arr, struct Points points, int N, mpfr_t nu, mpfr_t mu0,
             int (*index)(int)) {
  typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;
  mpreal *A_arr_mpreal;
  A_arr_mpreal = new mpreal[(points.boundary + points.interior)*N];

  // Fill A_arr with the coefficients for the matrix A
  generate_matrix(A_arr, points, N, nu, mu0, index);

  // Create an array of mpreals with the same values
  for (int i = 0; i < (points.boundary + points.interior)*N; i++)
    A_arr_mpreal[i] = mpreal(A_arr[i]);

  // Create the matrix A from the coefficients
  Map<MatrixXmp> A(A_arr_mpreal, points.boundary + points.interior, N);

  // Perform QR decomposition
  MatrixXmp Q(MatrixXmp::Identity(points.boundary + points.interior, N));
  HouseholderQR<MatrixXmp> qr(A);
  Q = qr.householderQ() * Q;

  // Find the smallest singular values
  BDCSVD<MatrixXmp> svd(Q.block(0, 0, points.interior, N));

  delete [] A_arr_mpreal;
  return svd.singularValues()(N-1);
}

void coefs_sigma(mpfr_t *coefs_arr, mpfr_t *A_arr, struct Points points,
                 int N, mpfr_t nu, mpfr_t mu0, int (*index)(int)) {
  typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;
  mpreal *A_arr_mpreal;
  A_arr_mpreal = new mpreal[(points.boundary + points.interior)*N];

  // Fill A_arr with the coefficients for the matrix A
  generate_matrix(A_arr, points, N, nu, mu0, index);

  // Create an array of mpreals with the same values
  for (int i = 0; i < (points.boundary + points.interior)*N; i++)
    A_arr_mpreal[i] = mpreal(A_arr[i]);

  // Create the matrix A from the coefficients
  Map<MatrixXmp> A(A_arr_mpreal, points.boundary + points.interior, N);

  // Perform QR decomposition
  MatrixXmp Q(MatrixXmp::Identity(points.boundary + points.interior, N));
  HouseholderQR<MatrixXmp> qr(A);
  Q = qr.householderQ() * Q;

  // Compute the right singular vector corresponding to the smallest
  // singular value
  BDCSVD<MatrixXmp> svd(Q.block(0, 0, points.interior, N), ComputeThinV);

  // Compute the coefficients
  MatrixXmp coefs = qr.solve(Q*(svd.matrixV().col(N-1)));;

  coefs /= coefs(0);

  for (int i = 0; i < N; i++)
    mpfr_set(coefs_arr[i], coefs(i).mpfr_srcptr(), MPFR_RNDN);

  delete [] A_arr_mpreal;
}
