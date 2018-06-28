#include "tools.h"
#include "generate-matrix.h"
#include "sigma.h"
#include "enclose.h"

#include "Eigen/Core"
#include "mpreal.h"
#include "time.h"

#include <iostream>
#include <iomanip>
#include <unistd.h>

using namespace std;
using namespace mpfr;

void minimize_sigma(mpfr_t *A_arr, struct Points points, mpfr_t *scaling, int N,
                    mpfr_t a, mpfr_t b, mpfr_t mu0, mpreal tol,
                    int (*index)(int)) {
  mpfr_t c, d;
  mpreal invphi, invphi2, h, yc, yd;
  int n;

  mpfr_init(c);
  mpfr_init(d);

  invphi = (sqrt(5) - 1)/2;
  invphi2 = (3 - sqrt(5))/2;

  h = mpreal(b) - mpreal(a);

  if (h <= tol)
    return;

  n = (int) ceil(log(tol/h)/log(invphi));
  mpfr_add(c, a, (invphi2*h).mpfr_srcptr(), MPFR_RNDN);
  mpfr_add(d, a, (invphi*h).mpfr_srcptr(), MPFR_RNDN);

  yc = sigma(A_arr, points, scaling, N, c, mu0, index);
  yd = sigma(A_arr, points, scaling, N, d, mu0, index);

  for (int k = 0; k < n; k++) {
    if (yc < yd) {
      mpfr_set(b, d, MPFR_RNDN);
      mpfr_set(d, c, MPFR_RNDN);
      yd = yc;
      h = invphi*h;
      mpfr_add(c, a, (invphi2*h).mpfr_srcptr(), MPFR_RNDN);
      yc = sigma(A_arr, points, scaling, N, c, mu0, index);
    } else {
      mpfr_set(a, c, MPFR_RNDN);
      mpfr_set(c, d, MPFR_RNDN);
      yc = yd;
      h = invphi*h;
      mpfr_add(d, a, (invphi*h).mpfr_srcptr(), MPFR_RNDN);
      yd = sigma(A_arr, points, scaling, N, d, mu0, index);
    }
  }

  if (yc < yd) {
    mpfr_set(b, d, MPFR_RNDN);
  } else {
    mpfr_set(a, c, MPFR_RNDN);
    mpfr_set(d, b, MPFR_RNDN);
  }

  mpfr_clear(c);
  mpfr_clear(d);
}

void particular_solution(struct Geometry geometry, int angles_coefs[],
                         mpfr_t *scaling, mpfr_t mu0, mpfr_t nu_guess,
                         int N, int (*index)(int), int output) {

  mpfr_t *A_arr, *coefs, *values;
  mpfr_t mu, nu_low, nu_upp, nu_step, nu, eps;
  mpreal s;
  int iterations;
  struct Points points, points_eigen;

  mpfr_init(mu);
  mpfr_init(nu_low);
  mpfr_init(nu_upp);
  mpfr_init(nu_step);
  mpfr_init(nu);
  mpfr_init(eps);

  points.boundary = 2*N;
  points.interior = 2*N;
  points_eigen.boundary = 500;
  points_eigen.interior = 0;

  points_init(points);
  points_init(points_eigen);

  A_arr = new mpfr_t[(points.boundary + points.interior)*N];
  coefs = new mpfr_t[N];
  values = new mpfr_t[points_eigen.boundary + points_eigen.interior];

  for (int i = 0; i < (points.boundary + points.interior)*N; i++) {
    mpfr_init(A_arr[i]);
  }
  for (int i = 0; i < N; i++) {
    mpfr_init(coefs[i]);
  }
  for (int i = 0; i < points_eigen.interior + points_eigen.boundary; i++) {
    mpfr_init(values[i]);
  }

  boundary(points, geometry);
  interior(points, geometry);

  if (output == 0) {
    // Plot the values of sigma
    iterations = 400;
    mpfr_sub_d(nu_low, nu_guess, 1e-4, MPFR_RNDN);
    mpfr_add_d(nu_upp, nu_guess, 1e-4, MPFR_RNDN);
    mpfr_sub(nu_step, nu_upp, nu_low, MPFR_RNDN);
    mpfr_div_si(nu_step, nu_step, iterations, MPFR_RNDN);

    cout << N << " " << iterations << endl;
    for (int i = 0; i < iterations; i++) {
      mpfr_set(nu, nu_step, MPFR_RNDN);
      mpfr_mul_si(nu, nu, i, MPFR_RNDN);
      mpfr_add(nu, nu, nu_low, MPFR_RNDN);
      s = sigma(A_arr, points, scaling, N, nu, mu0, index);
      cout << mpreal(nu) << " " << s << endl;
    }
  }

  if (output > 0) {
    // Find the value of nu that minimizes sigma
    mpfr_sub_d(nu_low, nu_guess, 1e-2, MPFR_RNDN);
    mpfr_add_d(nu_upp, nu_guess, 1e-2, MPFR_RNDN);
    minimize_sigma(A_arr, points, scaling, N, nu_low, nu_upp, mu0,
                   mpreal(1e-10), index);

    mpfr_add(nu, nu_low, nu_upp, MPFR_RNDN);
    mpfr_div_si(nu, nu, 2, MPFR_RNDN);

    cout << N << " " << mpreal(nu);
    if (output == 2)
      cout << " " << points_eigen.boundary << endl;
    else
      cout << endl;

    // Find the coefficients of the expansion
    coefs_sigma(coefs, A_arr, points, scaling, N, nu, mu0, index);

    if (output == 1)
      // Print the coefficients for the eigenfunction
      for (int i = 0; i < N; i++) {
        cout << mpreal(coefs[i]) << endl;
      }

    if (output == 2) {
      // Plotting the eigenfunction
      boundary(points_eigen, geometry);
      eigenfunction(values, coefs, points_eigen, N, nu, mu0, index);

      for (int i = 0; i < points_eigen.boundary; i++) {
        cout << mpreal(points_eigen.phis[i]) << " " << mpreal(values[i]) << endl;
      }
    }

    if (output == 3) {
      // Compute an enclosure for the eigenfunction This part is not
      // completely safe to rounding errors, the lower bound should be
      // rounded down and the upper bound up.
      enclose(eps, angles_coefs, coefs, N, nu, index);
      cout << mpreal(eps) << endl
           << mpreal(nu)*(1 + mpreal(nu))/(1 + mpreal(eps)) << " "
           << mpreal(nu)*(1 + mpreal(nu))/(1 - mpreal(eps)) << endl;
    }
  }
  /* Clear all variables */
  for (int i = 0; i < (points.boundary + points.interior)*N; i++) {
    mpfr_clear(A_arr[i]);
  }
  for (int i = 0; i < N; i++) {
    mpfr_clear(coefs[i]);
  }
  for (int i = 0; i < points_eigen.boundary + points_eigen.interior; i++) {
    mpfr_clear(values[i]);
  }
  delete [] A_arr;
  delete [] coefs;
  delete [] values;

  points_clear(points);
  points_clear(points_eigen);

  mpfr_clear(mu);
  mpfr_clear(nu_low);
  mpfr_clear(nu_upp);
  mpfr_clear(nu_step);
  mpfr_clear(nu);
  mpfr_clear(eps);
}

int index_function(int k) {
  return 2*k + 1;
}

int main(int argc, char *argv[]) {
  mpfr_t angles[3], *scaling;
  mpfr_t mu0, nu_guess, mu;
  int angles_coefs[6];
  int c, prec, output, N_beg, N_end, N_step;
  string usage;
  struct Geometry geometry;

  geometry_init(geometry);

  for (int i = 0; i < 3; i++) {
    mpfr_init(angles[i]);
  }
  mpfr_init(mu0);
  mpfr_init(nu_guess);
  mpfr_init(mu);

  srand(1);
  cout << setprecision(20);

  usage = "Usage: ./particular-solution [OPTION]... -- N1 D1 N2 D2 N3 D3\n\
Evaluate the method of particular solution for the spherical triangle given \n\
by the angles (pi*N1/D1, pi*N2/D2, pi*N3/D3). By default it uses N from 4 to 16\n\
with a step size of 2.\n\
Options are:\n\
  -n <value< - guess for the eigenvalue, used as midpoint (default 1.825757)\n\
  -p <value> - set the working precision to this (default 53)\n\
  -o <value> - set output type, valid values are 0, 1, 2.\n\
               0: Output data for plotting the values of sigma\n\
               1: Output the coefficients of the expansions\n\
               2: Output data for plotting the approximate eigenfunctions\n\
  -b <value> - start value for N (default 4)\n\
  -e <value> - end value for N (default 16)\n\
  -s <value> - step size for N (default 2)";

  // Set default values for parameters
  prec = 53;
  output = 2;
  N_beg = 4;
  N_end = 16;
  N_step = 2;
  mpfr_set_str(nu_guess, "1.825757", 10, MPFR_RNDN);

  while ((c = getopt (argc, argv, "n:p:o:b:e:s:")) != -1)
    switch(c) {
    case 'n':
      mpfr_set_str(nu_guess, optarg, 10, MPFR_RNDN);
      break;
    case 'p':
      prec = atoi(optarg);
      break;
    case 'o':
      output = atoi(optarg);
      break;
    case 'b':
      N_beg = atoi(optarg);
      break;
    case 'e':
      N_end = atoi(optarg);
      break;
    case 's':
      N_step = atoi(optarg);
      break;
    case '?':
      if (optopt == 'p')
        cerr <<"Option -" << char(optopt) << " requires an argument.\n\n"
             << endl;
      else if (isprint (optopt))
        cerr << "Unknown option `-" << char(optopt) << "'.\n\n" << endl;
      else
        cerr << "Unknown option character `" << char(optopt) << "'.\n\n"
             << endl;
      cerr << usage << endl;
      exit(0);
    default:
      abort ();
    }

  mpfr_set_default_prec(prec);

  mpfr_const_pi(angles[0], MPFR_RNDN);
  mpfr_const_pi(angles[1], MPFR_RNDN);
  mpfr_const_pi(angles[2], MPFR_RNDN);
  if (argc - optind > 1) {
    mpfr_mul_si(angles[0], angles[0], atoi(argv[optind]), MPFR_RNDN);
    mpfr_div_si(angles[0], angles[0], atoi(argv[optind + 1]), MPFR_RNDN);
    mpfr_set_si(mu0, -atoi(argv[optind + 1]), MPFR_RNDN);
    mpfr_div_si(mu0, mu0, atoi(argv[optind]), MPFR_RNDN);
    angles_coefs[0] = atoi(argv[optind]);
    angles_coefs[1] = atoi(argv[optind + 1]);
  } else {
    mpfr_mul_si(angles[0], angles[0], 2, MPFR_RNDN);
    mpfr_div_si(angles[0], angles[0], 3, MPFR_RNDN);
    mpfr_set_si(mu0, -3, MPFR_RNDN);
    mpfr_div_si(mu0, mu0, 2, MPFR_RNDN);
    angles_coefs[0] = 2;
    angles_coefs[1] = 3;
  }
  if (argc - optind > 3) {
    mpfr_mul_si(angles[1], angles[1], atoi(argv[optind + 2]), MPFR_RNDN);
    mpfr_div_si(angles[1], angles[1], atoi(argv[optind + 3]), MPFR_RNDN);
    angles_coefs[2] = atoi(argv[optind + 2]);
    angles_coefs[3] = atoi(argv[optind + 3]);
  } else {
    mpfr_mul_si(angles[1], angles[1], 2, MPFR_RNDN);
    mpfr_div_si(angles[1], angles[1], 3, MPFR_RNDN);
    angles_coefs[2] = 2;
    angles_coefs[3] = 3;
  }
  if (argc - optind > 5) {
    mpfr_mul_si(angles[2], angles[2], atoi(argv[optind + 4]), MPFR_RNDN);
    mpfr_div_si(angles[2], angles[2], atoi(argv[optind + 5]), MPFR_RNDN);
    angles_coefs[4] = atoi(argv[optind + 4]);
    angles_coefs[5] = atoi(argv[optind + 5]);
  } else {
    mpfr_mul_si(angles[2], angles[2], 2, MPFR_RNDN);
    mpfr_div_si(angles[2], angles[2], 3, MPFR_RNDN);
    angles_coefs[4] = 2;
    angles_coefs[5] = 3;
  }

  angles_to_vectors(geometry, angles);
  geometry.half_boundary = 1;

  scaling = new mpfr_t[N_end];
  for (int i = 0; i < N_end; i++) {
    mpfr_init(scaling[i]);
    mpfr_mul_si(mu, mu0, index_function(i), MPFR_RNDN);
    scale_norm(scaling[i], geometry.theta_bound, nu_guess, mu);
  }

  for (int N = N_beg; N <= N_end; N+=N_step) {
    particular_solution(geometry, angles_coefs, scaling, mu0, nu_guess, N,
                        index_function, output);
  }

  for (int i = 0; i < N_end; i++) {
    mpfr_clear(scaling[i]);
  }

  delete [] scaling;

  for (int i = 0; i < 3; i++) {
    mpfr_clear(angles[i]);
  }

  geometry_clear(geometry);

  mpfr_clear(mu0);
  mpfr_clear(nu_guess);
  mpfr_clear(mu);

  return 0;
}
