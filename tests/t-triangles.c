#include "particular_solution.h"
#include "arb.h"

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
  3.560131e-09,
  4.876337e-04,
  2.840797e-21,
  1.529307e-07,
  3.143688e-14,
  9.798216e-02,
  4.478798e-02,
  1.817899e-01,
  3.847273e-02
};

int main()
{
  arb_t nu_enclosure, ans, res;
  geom_t geometry;
  options_t options;
  slong prec;

  arb_init(nu_enclosure);
  arb_init(ans);
  arb_init(res);

  prec = 64;

  geom_init(geometry);

  for (int i = 0; i < NUM_TRIANGLES; i++)
  {
    options_default(options);
    get_domain(geometry, nu_enclosure, options, i);
    particular_solution_enclosure(nu_enclosure, geometry, options, prec);

    /* lambda = nu*(nu + 1) */
    arb_set(res, nu_enclosure);
    arb_add_si(ans, res, 1, 2*prec);
    arb_mul(res, res, ans, 2*prec);

    arb_set_str(ans, ans_str[i], 2*prec);

    flint_printf("Triangle %i\n", i);
    flint_printf("ans = "); arb_printn(ans, 25,  0); flint_printf("\n");
    flint_printf("res = "); arb_printn(res, 25,  0); flint_printf("\n");
    if (!arb_overlaps(ans, res))
    {
      flint_printf("FAILED TEST: Computed enclosure doesn't contain eigenvalue\n", i);
      flint_abort();
    }
    if (width[i] < mag_get_d(arb_radref(res)))
    {
      flint_printf("LOW PRECISION\n", i);
    }
    flint_printf("width goal = %e\n", width[i]);
    flint_printf("width res  = %e\n", mag_get_d(arb_radref(res)));
  }

  arb_clear(nu_enclosure);
  arb_clear(ans);
  arb_clear(res);

  geom_clear(geometry);

  flint_cleanup();

  return 0;
}
