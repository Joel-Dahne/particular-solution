#include "particular_solution.h"

#include <iostream>
#include <iomanip>
#include <unistd.h>

using namespace std;

int main(int argc, char *argv[]) {
  arb_t nu_enclosure;
  mpfr_t angles[3];
  mpfr_t nu_guess, nu_width, tol_rel, tol, tmp;
  double prec_factor;
  slong angles_coefs[6];
  int c, prec, output, N_beg, N_end, N_step;
  string usage;
  char nu_guess_str_default[] = "4.0631";
  char nu_width_str_default[] = "1e-2";
  char tol_rel_str_default[] = "1e-5";
  char *nu_guess_str, *nu_width_str, *tol_rel_str;
  geom_t geometry;
  int (*index_function)(int);
  particular_solution_opt_t options;

  particular_solution_opt_init(options);

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
               5: Output the validated enclosure of nu\n\
               6: Output the width of the validated enclosure of the eigenvalue\n\
               7: Output the width of the validated enclosure of nu\n\
  -b <value> - start value for N (default 4)\n\
  -e <value> - end value for N (default 16)\n\
  -s <value> - step size for N (default 2)\n\
  -h         - set this flag to use only half of the boundary\n\
  -i         - set this flag to use only odd indices in the expansion";

  // Set default values for parameters
  prec_factor = 1.2;
  prec = 53;
  output = 4;
  N_beg = 4;
  N_end = 16;
  N_step = 2;
  nu_guess_str = nu_guess_str_default;
  nu_width_str = nu_width_str_default;
  tol_rel_str = tol_rel_str_default;
  geometry->half_edge[0] = 0;
  geometry->half_edge[1] = 0;
  geometry->half_edge[2] = 0;
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
      geometry->half_edge[0] = 1;
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

  geom_init(geometry);

  for (int i = 0; i < 3; i++) {
    mpfr_init(angles[i]);
  }

  mpfr_init(nu_guess);
  mpfr_init(nu_width);
  mpfr_init(tol_rel);
  mpfr_init(tol);
  mpfr_init(tmp);

  arb_init(nu_enclosure);
  arb_set_str(nu_enclosure, nu_guess_str, prec);
  mag_set_d(arb_radref(nu_enclosure), atof(nu_width_str));

  mpfr_set_str(tol_rel, tol_rel_str, 10, MPFR_RNDN);

  for (int N = N_beg; N <= N_end; N+=N_step) {
    mpfr_mul(tol, tol, tol_rel, MPFR_RNDN);

    mpfr_log2(tmp, tol, MPFR_RNDN);
    mpfr_mul_d(tmp, tmp, -prec_factor, MPFR_RNDN);
    if (mpfr_get_si(tmp, MPFR_RNDN) > prec)
      prec = mpfr_get_si(tmp, MPFR_RNDN);

    mpfr_set_default_prec(prec);

    /* Recompute values with new precision */
    geom_set_angles(geometry, angles_coefs);

    particular_solution_enclosure(nu_enclosure, geometry, options, prec);
  }

  geom_clear(geometry);
  for (int i = 0; i < 3; i++) {
    mpfr_clear(angles[i]);
  }
  mpfr_clear(nu_guess);
  mpfr_init(nu_width);
  mpfr_clear(tol);
  mpfr_clear(tmp);

  return 0;
}
