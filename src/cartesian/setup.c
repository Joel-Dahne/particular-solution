#include "setup.h"

/* Get geometry and options for a number of default planar
 * triangles. */
void
get_domain(geom_t geometry, arb_t lambda_enclosure, options_t options,
           int triangle)
{
  slong angle[2];

  if (triangle == 0)
  {
    /* Equilateral triangle */

    /* Set up the coefficients for the angles */
    angle[0] = 1;
    angle[1] = 3;

    /* Set the side length */
    arb_set_si(geometry->side_length, 1);

    arb_set_d(lambda_enclosure, 52.637890139143245);
  }
  else if (triangle == 1)
  {
    /* Isosceles triangle */

    /* Set up the coefficients for the angles */
    angle[0] = 2;
    angle[1] = 3;

    /* Set the side length */
    arb_set_si(geometry->side_length, 1);

    arb_set_d(lambda_enclosure, 71.71);
  }
  else if (triangle == 2)
  {
    /* Isosceles triangle */

    /* Set up the coefficients for the angles */
    angle[0] = 2;
    angle[1] = 5;

    /* Set the side length */
    arb_set_si(geometry->side_length, 1);

    arb_set_d(lambda_enclosure, 48.598276257514954);
  }
  else if (triangle == 3)
  {
    /* Triangle with one angle 2pi/9 and the others about 3pi/9 and
     * 4pi/9 */

    /* Set up the coefficients for the angles */
    angle[0] = 2;
    angle[1] = 9;

    /* Set the side length */
    arb_set_d(geometry->side_length, 0.760995125438058);

    arb_set_d(lambda_enclosure, 101.84212553501129);
  }

  /* Set up the geometry */
  geom_set_angles(geometry, angle);

  /* These values are currently not used but needs to be set anyway */
  geometry->vertices[0] = 1;
  geometry->vertices[1] = 0;
  geometry->vertices[2] = 0;
  geometry->half_edge[0] = 0;
  geometry->half_edge[1] = 0;
  geometry->half_edge[2] = 0;


  /* Set starting enclosure of nu */
  mag_set_d(arb_radref(lambda_enclosure), 0.5);
}
