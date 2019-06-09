#include "setup.h"

void
particular_solution_opt_default(particular_solution_opt_t options)
{
  options->tol_relative = 1e-5;
  options->prec_factor = 1.2;
  options->N_beg = 4;
  options->N_end = 16;
  options->N_step = 2;
  options->plot_n = 0;
  options->output = 0;
  options->output_final = 0;
  options->output_time = 0;
}

/* Get geometry and options for a number of default triangles. */
void
get_triangle_defaults(geom_t geometry, arb_t nu_enclosure,
                      particular_solution_opt_t options, int triangle)
{
  slong angles[6];

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
    options->prec_factor = 2.0;

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
    options->prec_factor = 2.0;

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
    options->prec_factor = 2.0;

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
  mag_set_d(arb_radref(nu_enclosure), 0.5);
}
