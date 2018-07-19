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
#include <math.h>

using namespace std;
using namespace mpfr;

void minimize_sigma(mpfr_t nu, mpfr_t *A_arr, struct Points points, int N,
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

void particular_solution(struct Geometry geometry, int angles_coefs[],
                         mpfr_t mu0, mpfr_t nu_low, mpfr_t nu_upp, mpfr_t tol,
                         int N, int (*index)(int), int output) {

  mpfr_t *A_arr, *coefs, *values;
  mpfr_t mu, nu_step, nu, eps;
  mpreal s;
  int iterations;
  struct Points points, points_eigen;

  mpfr_init(mu);
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
    mpfr_sub(nu_step, nu_upp, nu_low, MPFR_RNDN);
    mpfr_div_si(nu_step, nu_step, iterations, MPFR_RNDN);
    cout << N << " " << iterations << endl;
    for (int i = 0; i < iterations; i++) {
      mpfr_set(nu, nu_step, MPFR_RNDN);
      mpfr_mul_si(nu, nu, i, MPFR_RNDN);
      mpfr_add(nu, nu, nu_low, MPFR_RNDN);
      s = sigma(A_arr, points, N, nu, mu0, index);
      cout << mpreal(nu) << " " << s << endl;
    }
  } else {
    // Find the value of nu that minimizes sigma
    minimize_sigma(nu, A_arr, points, N, nu_low, nu_upp, mu0, mpreal(tol),
                   index);

    cout << N << flush;
    if (output <= 3) {
      cout << " " << mpreal(nu);
      if (output == 3)
        cout << " " << points_eigen.boundary << endl;
      else
        cout << endl;
    }

    if (output >= 2) {
      // Find the coefficients of the expansion
      coefs_sigma(coefs, A_arr, points, N, nu, mu0, index);

      if (output == 2) {
        // Print the coefficients for the eigenfunction
        for (int i = 0; i < N; i++) {
          cout << mpreal(coefs[i]) << endl;
        }
      }

      if (output == 3) {
        // Plot the eigenfunction
        boundary(points_eigen, geometry);
        eigenfunction(values, coefs, points_eigen, N, nu, mu0, index);

        for (int i = 0; i < points_eigen.boundary; i++) {
          cout << mpreal(points_eigen.phis[i]) << " " << mpreal(values[i]) << endl;
        }
      }

      if (output == 4 || output == 5) {
        // Compute an enclosure of the eigenvalue. To be sure to get
        // correct output of the eigenvalue the output of it is
        // handled inside the function enclose.
        enclose(nu_low, nu_upp, angles_coefs, coefs, N, nu, index, output);
      }
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
  mpfr_clear(nu_step);
  mpfr_clear(nu);
  mpfr_clear(eps);
}

int index_function_odd(int k) {
  return 2*k + 1;
}

int index_function_all(int k) {
  return k + 1;
}

int main(int argc, char *argv[]) {
  mpfr_t angles[3];
  mpfr_t mu0, nu_guess, nu_width, nu_low, nu_upp, tol_rel, tol;
  double prec_factor;
  int angles_coefs[6];
  int c, prec, output, N_beg, N_end, N_step;
  string usage;
  char nu_guess_str_default[] = "4.0631";
  char nu_width_str_default[] = "1e-2";
  char tol_rel_str_default[] = "1e-5";
  char *nu_guess_str, *nu_width_str, *tol_rel_str;
  struct Geometry geometry;
  int (*index_function)(int);

  srand(1);

  usage = "Usage: ./particular-solution [OPTION]... -- N1 D1 N2 D2 N3 D3\n\
Evaluate the method of particular solution for the spherical triangle given \n\
by the angles (pi*N1/D1, pi*N2/D2, pi*N3/D3). By default it uses N from 4 to 16\n\
with a step size of 2.\n\
Options are:\n\
  -n <value> - guess for the eigenvalue, used as midpoint (default 4.0631)\n\
  -w <value> - width to use around the eigenvalue (default 1e-2)\n\
  -t <value> - relative, to the width, tolerance to use for minimization (default 1e-5)\n\
  -p <value> - set the minimum working precision to this (default 53)\n\
  -f <value> - set the factor used when computing the working precision needed\n\
               to achieve a given tolerance (default 1.2)\n\
  -o <value> - set output type, valid values are 0, 1, 2.\n\
               0: Output data for plotting the values of sigma\n\
               1: Output the nu value minimizing sigma\n\
               2: Output the coefficients of the expansions\n\
               3: Output data for plotting the approximate eigenfunctions\n\
               4: Output the validated enclosure of the eigenvalue\n\
               4: Output the validated enclosure of nu\n\
  -b <value> - start value for N (default 4)\n\
  -e <value> - end value for N (default 16)\n\
  -s <value> - step size for N (default 2)\n\
  -h         - set this flag to use only half of the boundary\n\
  -i         - set this flag to use only odd indices in the expansion";

  // Set default values for parameters
  prec_factor = 1.2;
  prec = 53;
  output = 1;
  N_beg = 4;
  N_end = 16;
  N_step = 2;
  nu_guess_str = nu_guess_str_default;
  nu_width_str = nu_width_str_default;
  tol_rel_str = tol_rel_str_default;
  geometry.half_boundary = 0;
  index_function = index_function_all;

  while ((c = getopt (argc, argv, "n:w:t:p:f:o:b:e:s:hi")) != -1)
    switch(c) {
    case 'n':
      nu_guess_str = optarg;
      break;
    case 'w':
      nu_width_str = optarg;
      break;
    case 't':
      tol_rel_str = optarg;
      break;
    case 'p':
      prec = atoi(optarg);
      break;
    case 'f':
      prec_factor = atof(optarg);
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
    case 'h':
      geometry.half_boundary = 1;
      break;
    case 'i':
      index_function = index_function_odd;
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

  if (argc - optind > 1) {
    angles_coefs[0] = atoi(argv[optind]);
    angles_coefs[1] = atoi(argv[optind + 1]);
  } else {
    angles_coefs[0] = 2;
    angles_coefs[1] = 3;
  }
  if (argc - optind > 3) {
    angles_coefs[2] = atoi(argv[optind + 2]);
    angles_coefs[3] = atoi(argv[optind + 3]);
  } else {
    angles_coefs[2] = 1;
    angles_coefs[3] = 4;
  }
  if (argc - optind > 5) {
    angles_coefs[4] = atoi(argv[optind + 4]);
    angles_coefs[5] = atoi(argv[optind + 5]);
  } else {
    angles_coefs[4] = 1;
    angles_coefs[5] = 2;
  }

  mpfr_set_default_prec(prec);
  cout << setprecision((int)ceil(prec*log10(2)));

  geometry_init(geometry);
  for (int i = 0; i < 3; i++) {
    mpfr_init(angles[i]);
  }
  mpfr_init(mu0);
  mpfr_init(nu_guess);
  mpfr_init(nu_width);
  mpfr_init(nu_low);
  mpfr_init(nu_upp);
  mpfr_init(tol_rel);
  mpfr_init(tol);

  mpfr_set_str(nu_guess, nu_guess_str, 10, MPFR_RNDN);
  mpfr_set_str(nu_width, nu_width_str, 10, MPFR_RNDN);
  mpfr_set_str(tol_rel, tol_rel_str, 10, MPFR_RNDN);

  mpfr_sub(nu_low, nu_guess, nu_width, MPFR_RNDN);
  mpfr_add(nu_upp, nu_guess, nu_width, MPFR_RNDN);

  for (int N = N_beg; N <= N_end; N+=N_step) {
    mpfr_sub(tol, nu_upp, nu_low, MPFR_RNDN);
    mpfr_mul(tol, tol, tol_rel, MPFR_RNDN);

    prec = mpfr_get_ui(max(-prec_factor*log2(mpreal(tol)), prec).mpfr_srcptr(),
                       MPFR_RNDN);
    mpfr_set_default_prec(prec);

    /* Round values to the new precision */
    mpfr_prec_round(nu_low, prec, MPFR_RNDN);
    mpfr_prec_round(nu_upp, prec, MPFR_RNDN);

    /* Recompute values with new precision */
    geometry_set_prec(geometry, prec);
    for (int i = 0; i < 3; i++) {
      mpfr_set_prec(angles[i], prec);
    }
    mpfr_set_prec(mu0, prec);

    mpfr_const_pi(angles[0], MPFR_RNDN);
    mpfr_const_pi(angles[1], MPFR_RNDN);
    mpfr_const_pi(angles[2], MPFR_RNDN);

    mpfr_mul_si(angles[0], angles[0], angles_coefs[0], MPFR_RNDN);
    mpfr_div_si(angles[0], angles[0], angles_coefs[1], MPFR_RNDN);
    mpfr_set_si(mu0, -angles_coefs[1], MPFR_RNDN);
    mpfr_div_si(mu0, mu0, angles_coefs[0], MPFR_RNDN);
    mpfr_mul_si(angles[1], angles[1], angles_coefs[2], MPFR_RNDN);
    mpfr_div_si(angles[1], angles[1], angles_coefs[3], MPFR_RNDN);
    mpfr_mul_si(angles[2], angles[2], angles_coefs[4], MPFR_RNDN);
    mpfr_div_si(angles[2], angles[2], angles_coefs[5], MPFR_RNDN);

    angles_to_vectors(geometry, angles);

    particular_solution(geometry, angles_coefs, mu0, nu_low, nu_upp, tol, N,
                        index_function, output);
  }

  geometry_clear(geometry);
  for (int i = 0; i < 3; i++) {
    mpfr_clear(angles[i]);
  }
  mpfr_clear(mu0);
  mpfr_clear(nu_guess);
  mpfr_init(nu_width);
  mpfr_init(nu_low);
  mpfr_init(nu_upp);
  mpfr_clear(tol);

  return 0;
}
