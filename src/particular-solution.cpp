#include "tools.h"
#include "generate-matrix.h"
#include "sigma.h"

#include "Eigen/Core"
#include "mpreal.h"
#include "time.h"

#include <iostream>
#include <iomanip>
#include <unistd.h>

using namespace std;
using namespace mpfr;

void minimize_sigma(mpfr_t *A_arr, mpfr_t *thetas, mpfr_t *phis,
                    mpfr_t *scaling, int boundary, int interior, int N,
                    mpfr_t a, mpfr_t b, mpfr_t mu0, mpreal tol,
                    int index_step) {
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

  yc = sigma(A_arr, thetas, phis, scaling, boundary, interior, N, c, mu0,
             index_step);
  yd = sigma(A_arr, thetas, phis, scaling, boundary, interior, N, d, mu0,
             index_step);

  for (int k = 0; k < n; k++) {
    if (yc < yd) {
      mpfr_set(b, d, MPFR_RNDN);
      mpfr_set(d, c, MPFR_RNDN);
      yd = yc;
      h = invphi*h;
      mpfr_add(c, a, (invphi2*h).mpfr_srcptr(), MPFR_RNDN);
      yc = sigma(A_arr, thetas, phis, scaling, boundary, interior, N, c, mu0,
                 index_step);
    } else {
      mpfr_set(a, c, MPFR_RNDN);
      mpfr_set(c, d, MPFR_RNDN);
      yc = yd;
      h = invphi*h;
      mpfr_add(d, a, (invphi*h).mpfr_srcptr(), MPFR_RNDN);
      yd = sigma(A_arr, thetas, phis, scaling, boundary, interior, N, d, mu0,
                 index_step);
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

void particular_solution(mpfr_t angles[], mpfr_t mu0, int N, int index_step,
                         int half_boundary, int output) {
  mpfr_t v1[3], v2[3], v3[3];
  mpfr_t *A_arr, *thetas, *phis, *scaling, *coefs, *thetas_eigen, *phis_eigen,
    *values;
  mpfr_t theta_bound, mu, nu_mid, nu_low, nu_upp, nu_step, nu;
  mpreal s;
  int num_boundary, num_interior, iterations, num_eigen;

  for (int i = 0; i < 3; i++) {
    mpfr_init(v1[i]);
    mpfr_init(v2[i]);
    mpfr_init(v3[i]);
  }

  mpfr_init(theta_bound);
  mpfr_init(mu);
  mpfr_init(nu_mid);
  mpfr_init(nu_low);
  mpfr_init(nu_upp);
  mpfr_init(nu_step);
  mpfr_init(nu);

  angles_to_vectors(v1, v2, v3, theta_bound, angles);

  num_boundary = 2*N;
  num_interior = 2*N;
  num_eigen = 500; //Should this depend on N?

  A_arr = new mpfr_t[(num_boundary + num_interior)*N];
  thetas = new mpfr_t[num_boundary + num_interior];
  phis = new mpfr_t[num_boundary + num_interior];
  scaling = new mpfr_t[N];
  coefs = new mpfr_t[N];
  thetas_eigen = new mpfr_t[num_eigen];
  phis_eigen = new mpfr_t[num_eigen];
  values = new mpfr_t[num_eigen];
  for (int i = 0; i < (num_boundary + num_interior)*N; i++) {
    mpfr_init(A_arr[i]);
  }
  for (int i = 0; i < num_boundary + num_interior; i++) {
    mpfr_init(thetas[i]);
    mpfr_init(phis[i]);
  }
  for (int i = 0; i < N; i++) {
    mpfr_init(scaling[i]);
    mpfr_init(coefs[i]);
  }
  for (int i = 0; i < num_eigen; i++) {
    mpfr_init(thetas_eigen[i]);
    mpfr_init(phis_eigen[i]);
    mpfr_init(values[i]);
  }

  boundary(thetas, phis, v2, v3, num_boundary, half_boundary);
  interior(thetas + num_boundary, phis + num_boundary, v1, v2, v3,
           num_interior);

  mpfr_set_d(nu_mid, 1.825757081, MPFR_RNDN);
  for (int i = 0; i < N; i++) {
    mpfr_mul_si(mu, mu0, index_step*i+1, MPFR_RNDN);
    scale_norm(scaling[i], theta_bound, nu_mid, mu);
  }

  if (output == 0) {
    // Plot the values of sigma
    iterations = 400;
    mpfr_sub_d(nu_low, nu_mid, 1e-4, MPFR_RNDN);
    mpfr_add_d(nu_upp, nu_mid, 1e-4, MPFR_RNDN);
    mpfr_sub(nu_step, nu_upp, nu_low, MPFR_RNDN);
    mpfr_div_si(nu_step, nu_step, iterations, MPFR_RNDN);

    cout << N << " " << iterations << endl;
    for (int i = 0; i < iterations; i++) {
      mpfr_set(nu, nu_step, MPFR_RNDN);
      mpfr_mul_si(nu, nu, i, MPFR_RNDN);
      mpfr_add(nu, nu, nu_low, MPFR_RNDN);
      s = sigma(A_arr, thetas, phis, scaling, num_boundary, num_interior, N, nu,
                mu0, index_step);
      cout << mpreal(nu) << " " << s << endl;
    }
  }

  if (output == 1 || output == 2) {
    // Find the value of nu that minimizes sigma
    mpfr_sub_d(nu_low, nu_mid, 1e-2, MPFR_RNDN);
    mpfr_add_d(nu_upp, nu_mid, 1e-2, MPFR_RNDN);
    minimize_sigma(A_arr, thetas, phis, scaling, num_boundary, num_interior, N,
                   nu_low, nu_upp, mu0, mpreal(1e-10), index_step);

    mpfr_add(nu, nu_low, nu_upp, MPFR_RNDN);
    mpfr_div_si(nu, nu, 2, MPFR_RNDN);

    cout << N << " " << mpreal(nu);
    if (output == 1)
      cout << " " << num_eigen << endl;
    else
      cout << endl;

    // Find the coefficients of the expansion
    coefs_sigma(coefs, A_arr, thetas, phis, scaling, num_boundary, num_interior,
                N, nu, mu0, index_step);

    if (output == 2)
      for (int i = 0; i < N; i++) {
        cout << mpreal(coefs[i]) << endl;
      }

    if (output == 1) {
      // Plotting the eigenfunction
      boundary(thetas_eigen, phis_eigen, v2, v3, num_eigen, half_boundary);
      eigenfunction(values, coefs, thetas_eigen, phis_eigen, num_eigen, N, nu,
                    mu0, index_step);
      for (int i = 0; i < num_eigen; i++) {
        cout << mpreal(phis_eigen[i]) << " " << mpreal(values[i]) << endl;
      }
    }
  }
  /* Clear all variables */
  for (int i = 0; i < (num_boundary + num_interior)*N; i++) {
    mpfr_clear(A_arr[i]);
  }
  for (int i = 0; i < num_boundary + num_interior; i++) {
    mpfr_clear(thetas[i]);
    mpfr_clear(phis[i]);
  }
  for (int i = 0; i < N; i++) {
    mpfr_clear(scaling[i]);
    mpfr_clear(coefs[i]);
  }
  for (int i = 0; i < num_eigen; i++) {
    mpfr_clear(thetas_eigen[i]);
    mpfr_clear(phis_eigen[i]);
    mpfr_clear(values[i]);
  }
  delete [] A_arr;
  delete [] thetas;
  delete [] phis;
  delete [] scaling;
  delete [] coefs;
  delete [] thetas_eigen;
  delete [] phis_eigen;
  delete [] values;

  for (int i = 0; i < 3; i++) {
    mpfr_clear(v1[i]);
    mpfr_clear(v2[i]);
    mpfr_clear(v3[i]);
  }
  mpfr_clear(theta_bound);
  mpfr_clear(mu);
  mpfr_clear(nu_mid);
  mpfr_clear(nu_low);
  mpfr_clear(nu_upp);
  mpfr_clear(nu_step);
  mpfr_clear(nu);
}

int main(int argc, char *argv[]) {
  mpfr_t angles[3];
  mpfr_t mu0;
  int c, prec, output, N_beg, N_end, N_step, index_step, half_boundary;
  string usage;

  for (int i = 0; i < 3; i++) {
    mpfr_init(angles[i]);
  }
  mpfr_init(mu0);

  srand(1);
  cout << setprecision(20);

  usage = "Usage: ./particular-solution [OPTION]... -- N_begin N_end\n\
Evaluate the method of particular solution with the number of therms in\n\
the expansion varying from N_begin to N_end with a default step size of 2\n\
Options are:\n\
  -p <value> - set the working precision to this (default 53)\n\
  -o <value> - set output type, valid values are 0, 1, 2.\n\
               0: Output data for plotting the values of sigma\n\
               1: Output data for plotting the approximate eigenfunctions\n\
               2: Output the coefficients of the expansions\n\
  -s <value> - step size for N (default 2)";

  // Set default values for parameters
  prec = 53;
  output = 1;
  N_beg = 4;
  N_end = 16;
  N_step = 2;

  while ((c = getopt (argc, argv, "p:o:s:")) != -1)
    switch(c) {
    case 'p':
      prec = atoi(optarg);
      break;
    case 'o':
      output = atoi(optarg);
      break;
    case 's':
      N_step = atoi(optarg);
      break;
    case '?':
      if (optopt == 'p')
        cerr <<"Option -" << char(optopt) << " requires an argument.\n\n" << endl;
      else if (isprint (optopt))
        cerr << "Unknown option `-" << char(optopt) << "'.\n\n" << endl;
      else
        cerr << "Unknown option character `" << char(optopt) << "'.\n\n" << endl;
      cerr << usage << endl;
      exit(0);
    default:
      abort ();
    }

  if (argc - optind > 0)
    N_beg = atoi(argv[optind]);
  if (argc - optind > 1)
    N_end = atoi(argv[optind + 1]);

  mpfr_set_default_prec(prec);

  index_step = 1;
  half_boundary = 0;

  mpfr_const_pi(angles[0], MPFR_RNDN);
  mpfr_mul_si(angles[0], angles[0], 2, MPFR_RNDN);
  mpfr_div_si(angles[0], angles[0], 3, MPFR_RNDN);
  mpfr_const_pi(angles[1], MPFR_RNDN);
  mpfr_mul_si(angles[1], angles[1], 2, MPFR_RNDN);
  mpfr_div_si(angles[1], angles[1], 3, MPFR_RNDN);
  mpfr_const_pi(angles[2], MPFR_RNDN);
  mpfr_mul_si(angles[2], angles[2], 2, MPFR_RNDN);
  mpfr_div_si(angles[2], angles[2], 3, MPFR_RNDN);

  mpfr_set_str(mu0, "-1.5", 10, MPFR_RNDN);

  for (int N = N_beg; N <= N_end; N+=N_step) {
    particular_solution(angles, mu0, N, index_step, half_boundary, output);
  }

  for (int i = 0; i < 3; i++) {
    mpfr_clear(angles[i]);
  }
  mpfr_clear(mu0);

  return 0;
}
