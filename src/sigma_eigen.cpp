#include "unsupported/Eigen/MPRealSupport"
#include "Eigen/QR"
#include "Eigen/SVD"
#include "mpreal.h"

using namespace mpfr;
using namespace Eigen;

typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;

void sigma_eigen(mpfr_t res, mpfr_t *A_mpfr, int interior, int rows, int N)
{
  mpreal *A_mpreal;
  MatrixXmp A, Q;
  HouseholderQR<MatrixXmp> qr;
  BDCSVD<MatrixXmp> svd;

  A_mpreal = new mpreal[rows*N];
  for (int i = 0; i < rows*N; i++)
    A_mpreal[i] = mpreal(A_mpfr[i]);

  // Create the matrix A from the coefficients
  A = Map<MatrixXmp>(A_mpreal, rows, N);

  // Perform QR decomposition
  Q = MatrixXmp::Identity(rows, N);
  qr = HouseholderQR<MatrixXmp>(A);
  Q = qr.householderQ() * Q;

  // Find the smallest singular value
  svd = BDCSVD<MatrixXmp>(Q.block(0, 0, interior, N));
  mpfr_set(res, svd.singularValues()(N-1).mpfr_srcptr(), MPFR_RNDN);

  delete [] A_mpreal;
}

void coefs_sigma_eigen(mpfr_t *coefs_mpfr, mpfr_t *A_mpfr, int interior,
                       int rows, int N)
{
  mpreal *A_mpreal;
  MatrixXmp A, Q, coefs;
  HouseholderQR<MatrixXmp> qr;
  BDCSVD<MatrixXmp> svd;

  A_mpreal = new mpreal[rows*N];
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
  svd = BDCSVD<MatrixXmp>(Q.block(0, 0, interior, N), ComputeThinV);

  // Compute the coefficients
  coefs = qr.solve(Q*(svd.matrixV().col(N-1)));;

  coefs /= coefs(0);

  for (int i = 0; i < N; i++)
    mpfr_set(coefs_mpfr[i], coefs(i).mpfr_srcptr(), MPFR_RNDN);
}
