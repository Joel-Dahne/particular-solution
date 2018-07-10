#include "tools.h"

#include "mpfr.h"
#include "stdlib.h"

void
points_init(struct Points & points)
{
  points.thetas = new mpfr_t[(points.boundary + points.interior)];
  points.phis = new mpfr_t[(points.boundary + points.interior)];
  for (int i = 0; i < points.interior + points.boundary; i++) {
    mpfr_init(points.thetas[i]);
    mpfr_init(points.phis[i]);
  }
}

void
points_clear(struct Points & points)
{
  for (int i = 0; i < points.interior + points.boundary; i++) {
    mpfr_clear(points.thetas[i]);
    mpfr_clear(points.phis[i]);
  }
  delete [] points.thetas;
  delete [] points.phis;
}

void
geometry_init(struct Geometry & geometry)
{
  for (int i = 0; i < 3; i++) {
    mpfr_init(geometry.v1[i]);
    mpfr_init(geometry.v2[i]);
    mpfr_init(geometry.v3[i]);
  }
  mpfr_init(geometry.theta_bound);
}

void
geometry_clear(struct Geometry & geometry)
{
  for (int i = 0; i < 3; i++) {
    mpfr_clear(geometry.v1[i]);
    mpfr_clear(geometry.v2[i]);
    mpfr_clear(geometry.v3[i]);
  }
  mpfr_clear(geometry.theta_bound);
}

void
angles_to_vectors(struct Geometry &geometry, mpfr_t *angles)
{
  mpfr_t S, tmp1, tmp2;

  mpfr_init(S);
  mpfr_init(tmp1);
  mpfr_init(tmp2);

  /* Compute the vectors for the spherical triangle */
  /* S = (angles[0] + angles[1] + angles[2])/2 */
  mpfr_add(S, angles[0], angles[1], MPFR_RNDN);
  mpfr_add(S, S, angles[2], MPFR_RNDN);
  mpfr_div_si(S, S, 2, MPFR_RNDN);

  /* v1 = [0, 0, 1] */
  mpfr_set_si(geometry.v1[0], 0, MPFR_RNDN);
  mpfr_set_si(geometry.v1[1], 0, MPFR_RNDN);
  mpfr_set_si(geometry.v1[2], 1, MPFR_RNDN);

  /* Compute theta for v2*/
  mpfr_cos(tmp1, S, MPFR_RNDN);

  mpfr_sub(tmp2, S, angles[1], MPFR_RNDN);
  mpfr_cos(tmp2, tmp2, MPFR_RNDN);
  mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);

  mpfr_sin(tmp2, angles[0], MPFR_RNDN);
  mpfr_div(tmp1, tmp1, tmp2, MPFR_RNDN);

  mpfr_sin(tmp2, angles[2], MPFR_RNDN);
  mpfr_div(tmp1, tmp1, tmp2, MPFR_RNDN);

  mpfr_neg(tmp1, tmp1, MPFR_RNDN);
  mpfr_sqrt(tmp1, tmp1, MPFR_RNDN);
  mpfr_asin(tmp1, tmp1, MPFR_RNDN);
  mpfr_mul_si(tmp1, tmp1, 2, MPFR_RNDN);

  /* Add theta to theta_bound */
  mpfr_set(geometry.theta_bound, tmp1, MPFR_RNDN);

  /* Compute v2 from knowing theta */
  mpfr_sin(geometry.v2[0], tmp1, MPFR_RNDN);
  mpfr_set_si(geometry.v2[1], 0, MPFR_RNDN);
  mpfr_cos(geometry.v2[2], tmp1, MPFR_RNDN);

  /* Compute theta for v3*/
  mpfr_cos(tmp1, S, MPFR_RNDN);

  mpfr_sub(tmp2, S, angles[2], MPFR_RNDN);
  mpfr_cos(tmp2, tmp2, MPFR_RNDN);
  mpfr_mul(tmp1, tmp1, tmp2, MPFR_RNDN);

  mpfr_sin(tmp2, angles[0], MPFR_RNDN);
  mpfr_div(tmp1, tmp1, tmp2, MPFR_RNDN);

  mpfr_sin(tmp2, angles[1], MPFR_RNDN);
  mpfr_div(tmp1, tmp1, tmp2, MPFR_RNDN);

  mpfr_neg(tmp1, tmp1, MPFR_RNDN);
  mpfr_sqrt(tmp1, tmp1, MPFR_RNDN);
  mpfr_asin(tmp1, tmp1, MPFR_RNDN);
  mpfr_mul_si(tmp1, tmp1, 2, MPFR_RNDN);

  /* Add theta to theta_bound and divide by 2*/
  mpfr_add(geometry.theta_bound, geometry.theta_bound, tmp1, MPFR_RNDN);
  mpfr_div_si(geometry.theta_bound, geometry.theta_bound, 2, MPFR_RNDN);

  /* Compute v3 from knowing theta */
  mpfr_sin(tmp2, tmp1, MPFR_RNDN);
  mpfr_cos(geometry.v3[0], angles[0], MPFR_RNDN);
  mpfr_mul(geometry.v3[0], geometry.v3[0], tmp2, MPFR_RNDN);
  mpfr_sin(geometry.v3[1], angles[0], MPFR_RNDN);
  mpfr_mul(geometry.v3[1], geometry.v3[1], tmp2, MPFR_RNDN);
  mpfr_cos(geometry.v3[2], tmp1, MPFR_RNDN);

  mpfr_clear(S);
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);
}

void
boundary(struct Points points, struct Geometry g)
{
  mpfr_t arcx, arcy, arcz, t, norm, tmp;
  int n;

  mpfr_init(arcx);
  mpfr_init(arcy);
  mpfr_init(arcz);
  mpfr_init(t);
  mpfr_init(norm);
  mpfr_init(tmp);

  n = points.boundary;

  for (int i = 0; i < n; i++) {
    mpfr_set_si(t, i + 1, MPFR_RNDN);
    if (g.half_boundary)
      mpfr_div_si(t, t, 2*n, MPFR_RNDN);
    else
      mpfr_div_si(t, t, n + 1, MPFR_RNDN);

    mpfr_sub(arcx, g.v3[0], g.v2[0], MPFR_RNDN);
    mpfr_mul(arcx, arcx, t, MPFR_RNDN);
    mpfr_add(arcx, arcx, g.v2[0], MPFR_RNDN);

    mpfr_sub(arcy, g.v3[1], g.v2[1], MPFR_RNDN);
    mpfr_mul(arcy, arcy, t, MPFR_RNDN);
    mpfr_add(arcy, arcy, g.v2[1], MPFR_RNDN);

    mpfr_sub(arcz, g.v3[2], g.v2[2], MPFR_RNDN);
    mpfr_mul(arcz, arcz, t, MPFR_RNDN);
    mpfr_add(arcz, arcz, g.v2[2], MPFR_RNDN);

    mpfr_sqr(norm, arcx, MPFR_RNDN);

    mpfr_sqr(tmp, arcy, MPFR_RNDN);
    mpfr_add(norm, norm, tmp, MPFR_RNDN);

    mpfr_sqr(tmp, arcz, MPFR_RNDN);
    mpfr_add(norm, norm, tmp, MPFR_RNDN);

    mpfr_sqrt(norm, norm, MPFR_RNDN);

    mpfr_div(arcx, arcx, norm, MPFR_RNDN);
    mpfr_div(arcy, arcy, norm, MPFR_RNDN);
    mpfr_div(arcz, arcz, norm, MPFR_RNDN);

    mpfr_acos(points.thetas[i], arcz, MPFR_RNDN);
    mpfr_atan2(points.phis[i], arcy, arcx, MPFR_RNDN);
  }

  mpfr_clear(arcx);
  mpfr_clear(arcy);
  mpfr_clear(arcz);
  mpfr_clear(t);
  mpfr_clear(norm);
  mpfr_clear(tmp);
}

void
interior(struct Points points, struct Geometry g)
{
  mpfr_t s, t, x, y, z, norm, tmp;

  mpfr_init(x);
  mpfr_init(y);
  mpfr_init(z);
  mpfr_init(s);
  mpfr_init(t);
  mpfr_init(norm);
  mpfr_init(tmp);

  for (int i = points.boundary; i < points.boundary + points.interior; i++)
    {
        /* We take s and t random with s in [0, 1) and s < 1 - t*/
      mpfr_set_d(s, ((double)rand()/RAND_MAX), MPFR_RNDN);
      mpfr_set_d(t, ((double)rand()/RAND_MAX), MPFR_RNDN);

        /* Set tmp to 1 - t */
        mpfr_sub_si(tmp, t, 1, MPFR_RNDN);
        mpfr_neg(tmp, tmp, MPFR_RNDN);
        /* If s < 1 - t set s = 1 - t, t = 1 - s to make them satisfy
           s < 1 - t. */
        if (mpfr_cmp(s, tmp) >= 0)
        {
            /* Set t to 1 - s */
            mpfr_sub_si(t, s, 1, MPFR_RNDN);
            mpfr_neg(t, t, MPFR_RNDN);

            /* Set s to 1 - t */
            mpfr_set(s, tmp, MPFR_RNDN);
        }

        mpfr_sub(x, g.v2[0], g.v1[0], MPFR_RNDN);
        mpfr_mul(x, x, s, MPFR_RNDN);
        mpfr_sub(tmp, g.v3[0], g.v1[0], MPFR_RNDN);
        mpfr_mul(tmp, tmp, t, MPFR_RNDN);
        mpfr_add(x, x, tmp, MPFR_RNDN);
        mpfr_add(x, x, g.v1[0], MPFR_RNDN);

        mpfr_sub(y, g.v2[1], g.v1[1], MPFR_RNDN);
        mpfr_mul(y, y, s, MPFR_RNDN);
        mpfr_sub(tmp, g.v3[1], g.v1[1], MPFR_RNDN);
        mpfr_mul(tmp, tmp, t, MPFR_RNDN);
        mpfr_add(y, y, tmp, MPFR_RNDN);
        mpfr_add(y, y, g.v1[1], MPFR_RNDN);

        mpfr_sub(z, g.v2[2], g.v1[2], MPFR_RNDN);
        mpfr_mul(z, z, s, MPFR_RNDN);
        mpfr_sub(tmp, g.v3[2], g.v1[2], MPFR_RNDN);
        mpfr_mul(tmp, tmp, t, MPFR_RNDN);
        mpfr_add(z, z, tmp, MPFR_RNDN);
        mpfr_add(z, z, g.v1[2], MPFR_RNDN);

        mpfr_sqr(norm, x, MPFR_RNDN);

        mpfr_sqr(tmp, y, MPFR_RNDN);
        mpfr_add(norm, norm, tmp, MPFR_RNDN);

        mpfr_sqr(tmp, z, MPFR_RNDN);
        mpfr_add(norm, norm, tmp, MPFR_RNDN);

        mpfr_sqrt(norm, norm, MPFR_RNDN);

        mpfr_div(x, x, norm, MPFR_RNDN);
        mpfr_div(y, y, norm, MPFR_RNDN);
        mpfr_div(z, z, norm, MPFR_RNDN);

        mpfr_acos(points.thetas[i], z, MPFR_RNDN);
        mpfr_atan2(points.phis[i], y, x, MPFR_RNDN);
    }

  mpfr_clear(x);
  mpfr_clear(y);
  mpfr_clear(z);
  mpfr_clear(s);
  mpfr_clear(t);
  mpfr_clear(norm);
  mpfr_clear(tmp);
}
