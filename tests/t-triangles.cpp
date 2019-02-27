#include "particular_solution.h"
#include "arb.h"

/* ------------------------------------------------------------------------- */
/*  Different triangles                                                      */
/* ------------------------------------------------------------------------- */

void
get_triangle(geom_t geometry, slong angles[], arb_t nu_enclosure,
             particular_solution_opt_t options, int triangle)
{
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
    arb_set_d(options->prec_factor, 2);

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 0;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 3.056691018);
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
    arb_set_d(options->prec_factor, 2);

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 0;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;


    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 3.240902298);
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
    arb_set_d(options->prec_factor, 2);

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 0;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 4.063109028);
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

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 1;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;


    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 4.143210850);
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

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 0;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 4.470604591);
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
    arb_set_d(options->prec_factor, 2);

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 1;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 6.525663100);
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

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 1;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 1.624084509);
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

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 1;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 1.825757081);
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

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 0;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 2.047890892);
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

    /* Set which edges are to be used and for which we only use half
     * of the boundary */
    geometry->vertices[0] = 1;
    geometry->vertices[1] = 0;
    geometry->vertices[2] = 0;
    geometry->half_edge[0] = 1;
    geometry->half_edge[1] = 0;
    geometry->half_edge[2] = 0;

    /* Set midpoint for nu */
    arb_set_d(nu_enclosure, 2.150869291);
  }

  /* Set up the geometry */
  geom_set_angles(geometry, angles);

  /* Set starting enclosure of nu */
  mag_set_d(arb_radref(nu_enclosure), 1e-2);
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

const double width[NUM_TRIANGLES] =
{
  1.850986e-05,
  3.518423e-09,
  4.830439e-04,
  3.534395e-18,
  1.529171e-07,
  3.143821e-14,
  4.258170e-02,
  4.445406e-02,
  5.105783e-02,
  3.847281e-02
};

int main()
{
  arb_t nu_enclosure, ans, res;
  geom_t geometry;
  particular_solution_opt_t options;
  slong angles[6];
  slong prec;

  arb_init(nu_enclosure);
  arb_init(ans);
  arb_init(res);

  prec = 64;

  geom_init(geometry);

  particular_solution_opt_init(options);

  for (int i = 0; i < NUM_TRIANGLES; i++)
  {
    particular_solution_opt_default(options);
    get_triangle(geometry, angles, nu_enclosure, options, i);
    particular_solution_enclosure(nu_enclosure, geometry, options, prec);

    arb_set(res, nu_enclosure);
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
    else if (width[i] < mag_get_d(arb_radref(res)))
    {
      flint_printf("LOW PRECISION (Triangle %i)\n", i);
      flint_printf("ans = %e\n", width[i]);
      flint_printf("res = %e\n", mag_get_d(arb_radref(res)));
    }

    flint_printf("Triangle %i\n", i);
    flint_printf("ans = "); arb_printn(ans, 20,  0); flint_printf("\n");
    flint_printf("res = "); arb_printn(res, 20,  0); flint_printf("\n");
    flint_printf("width goal = %e\n", width[i]);
    flint_printf("width res  = %e\n", mag_get_d(arb_radref(res)));
  }

  arb_clear(nu_enclosure);
  arb_clear(ans);
  arb_clear(res);

  geom_clear(geometry);

  particular_solution_opt_clear(options);

  flint_cleanup();

  return 0;
}
