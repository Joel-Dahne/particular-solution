#include "particular_solution.h"
#include "arb.h"

/* ------------------------------------------------------------------------- */
/*  Different triangles                                                      */
/* ------------------------------------------------------------------------- */

void
get_triangle(geom_t geometry, int angles[], mpfr_t nu_low, mpfr_t nu_upp,
             particular_solution_opt_t options, int triangle, int prec)
{
  /* Set default precision */
  mpfr_set_default_prec(prec);

  if (triangle == 0)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 3;
    angles[1] = 4;
    angles[2] = 1;
    angles[3] = 3;
    angles[4] = 1;
    angles[5] = 2;

    /* Set custom options */
    options->prec_factor = 2;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 0;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 3.056691018, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 1)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 2;
    angles[1] = 3;
    angles[2] = 1;
    angles[3] = 3;
    angles[4] = 1;
    angles[5] = 2;

    /* Set custom options */
    options->prec_factor = 2;
    options->index_function = index_function_all;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 0;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 3.240902298, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 2)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 2;
    angles[1] = 3;
    angles[2] = 1;
    angles[3] = 4;
    angles[4] = 1;
    angles[5] = 2;

    /* Set custom options */
    options->prec_factor = 2;
    options->index_function = index_function_all;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 0;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 4.063109028, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 3)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 2;
    angles[1] = 3;
    angles[2] = 1;
    angles[3] = 3;
    angles[4] = 1;
    angles[5] = 3;

    /* Set custom options */
    options->prec_factor = 1.2;
    options->index_function = index_function_all;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 1;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 4.143210850, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 4)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 3;
    angles[1] = 4;
    angles[2] = 1;
    angles[3] = 4;
    angles[4] = 1;
    angles[5] = 3;

    /* Set custom options */
    options->prec_factor = 1.2;
    options->index_function = index_function_all;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 0;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 4.470604591, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 5)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 2;
    angles[1] = 3;
    angles[2] = 1;
    angles[3] = 4;
    angles[4] = 1;
    angles[5] = 4;

    /* Set custom options */
    options->prec_factor = 1.2;
    options->index_function = index_function_odd;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 1;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 6.525663100, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 6)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 2;
    angles[1] = 3;
    angles[2] = 3;
    angles[3] = 4;
    angles[4] = 3;
    angles[5] = 4;

    /* Set custom options */
    options->prec_factor = 1.2;
    options->index_function = index_function_odd;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 1;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 1.624084509, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 7)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 2;
    angles[1] = 3;
    angles[2] = 2;
    angles[3] = 3;
    angles[4] = 2;
    angles[5] = 3;

    /* Set custom options */
    options->prec_factor = 1.2;
    options->index_function = index_function_odd;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 1;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 1.825757081, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 8)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 1;
    angles[1] = 2;
    angles[2] = 2;
    angles[3] = 3;
    angles[4] = 3;
    angles[5] = 4;

    /* Set custom options */
    options->prec_factor = 1.2;
    options->index_function = index_function_all;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 0;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 2.047890892, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }
  else if (triangle == 9)
  {
    /* Set up the coefficients for the angles */
    angles[0] = 1;
    angles[1] = 2;
    angles[2] = 2;
    angles[3] = 3;
    angles[4] = 2;
    angles[5] = 3;

    /* Set custom options */
    options->prec_factor = 1.2;
    options->index_function = index_function_odd;

    /* Set if only half of the boundary is to be used or not */
    geometry->half_boundary = 1;

    /* Set midpoint for nu */
    mpfr_set_prec(nu_low, prec);
    mpfr_set_prec(nu_upp, prec);
    mpfr_set_d(nu_low, 2.150869291, MPFR_RNDN);
    mpfr_set(nu_upp, nu_low, MPFR_RNDN);
  }

  /* Set up the geometry */
  geom_set(geometry, angles, prec);

  /* Set starting enclosure of nu */
  mpfr_sub_d(nu_low, nu_low, 1e-2, MPFR_RNDN);
  mpfr_add_d(nu_upp, nu_upp, 1e-2, MPFR_RNDN);
}

#define NUM_TRIANGLES 10

const char * ans_str[NUM_TRIANGLES] =
{
  "12.40005165284337790528605341 +/- 4.96e-27",
  "13.7443552132132318354011215921380207828066502596318748941363320689579830254389619211598920192317 +/- 7.13e-95",
  "20.5719735379847305566258421533 +/- 6.17e-29",
  "21.309407630190445258953481441230517778336842577146716613113142418206238547040233941912302059567611577883829836706377598939726916941225413300936673580274916786587 +/- 2.21e-160",
  "24.4569137962991116944804381447726828996080 +/- 6.73e-41",
  "49.109945263284609919670343151508268353698425615333956068479546500637275248339988486176558994445206617439284515387218370698834970763 +/- 3.91e-130",
  "4.3 +/- 0.0813",
  "5.16 +/- 7.74e-3",
  "6.2 +/- 0.0849",
  "6.78 +/- 8.61e-3"
};

int main()
{
  arb_t ans, res;
  mpfr_t nu_low, nu_upp;
  int angles[6];
  geom_t geometry;
  particular_solution_opt_t options;
  int prec;

  arb_init(ans);
  arb_init(res);

  prec = 64;

  mpfr_init(nu_low);
  mpfr_init(nu_upp);

  geom_init(geometry);

  particular_solution_opt_init(options);

  for (int i = 0; i < NUM_TRIANGLES; i++)
  {
    get_triangle(geometry, angles, nu_low, nu_upp, options, i, prec);
    particular_solution_enclosure(geometry, angles, nu_low, nu_upp,
                                  options);
    arb_set_interval_mpfr(res, nu_low, nu_upp, prec);
    arb_add_si(ans, res, 1, prec);
    arb_mul(res, res, ans, prec);

    arb_set_str(ans, ans_str[i], prec);

    if (!arb_overlaps(ans, res))
    {
      flint_printf("FAIL (Triangle %i)\n", i);
      flint_printf("ans = "); arb_printn(ans, 20,  0); flint_printf("\n");
      flint_printf("res = "); arb_printn(res, 20,  0); flint_printf("\n");
      flint_abort();
    }

    flint_printf("ans = "); arb_printn(ans, 20,  0); flint_printf("\n");
    flint_printf("res = "); arb_printn(res, 20,  0); flint_printf("\n");
  }

  arb_clear(ans);
  arb_clear(res);

  mpfr_clear(nu_low);
  mpfr_clear(nu_upp);

  geom_clear(geometry);

  flint_cleanup();

  return 0;
}
