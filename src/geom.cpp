#include "geom.h"

void
geom_init(geom_t g)
{
  g->v1 = new mpfr_t[3];
  g->v2 = new mpfr_t[3];
  g->v3 = new mpfr_t[3];
  for (int i = 0; i < 3; i++) {
    mpfr_init((g->v1)[i]);
    mpfr_init(g->v2[i]);
    mpfr_init(g->v3[i]);
  }
  mpfr_init(g->theta_bound);
}

void
geom_clear(geom_t g)
{
  for (int i = 0; i < 3; i++) {
    mpfr_clear(g->v1[i]);
    mpfr_clear(g->v2[i]);
    mpfr_clear(g->v3[i]);
  }
  mpfr_clear(g->theta_bound);
  delete [] g->v1;
  delete [] g->v2;
  delete [] g->v3;
}

void
geom_set_prec(geom_t g, mpfr_prec_t prec)
{
  for (int i = 0; i < 3; i++) {
    mpfr_set_prec(g->v1[i], prec);
    mpfr_set_prec(g->v2[i], prec);
    mpfr_set_prec(g->v3[i], prec);
  }
  mpfr_set_prec(g->theta_bound, prec);
}

void
geom_set(geom_t g, mpfr_t *angles)
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
  mpfr_set_si(g->v1[0], 0, MPFR_RNDN);
  mpfr_set_si(g->v1[1], 0, MPFR_RNDN);
  mpfr_set_si(g->v1[2], 1, MPFR_RNDN);

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
  mpfr_set(g->theta_bound, tmp1, MPFR_RNDN);

  /* Compute v2 from knowing theta */
  mpfr_sin(g->v2[0], tmp1, MPFR_RNDN);
  mpfr_set_si(g->v2[1], 0, MPFR_RNDN);
  mpfr_cos(g->v2[2], tmp1, MPFR_RNDN);

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
  mpfr_add(g->theta_bound, g->theta_bound, tmp1, MPFR_RNDN);
  mpfr_div_si(g->theta_bound, g->theta_bound, 2, MPFR_RNDN);

  /* Compute v3 from knowing theta */
  mpfr_sin(tmp2, tmp1, MPFR_RNDN);
  mpfr_cos(g->v3[0], angles[0], MPFR_RNDN);
  mpfr_mul(g->v3[0], g->v3[0], tmp2, MPFR_RNDN);
  mpfr_sin(g->v3[1], angles[0], MPFR_RNDN);
  mpfr_mul(g->v3[1], g->v3[1], tmp2, MPFR_RNDN);
  mpfr_cos(g->v3[2], tmp1, MPFR_RNDN);

  mpfr_clear(S);
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);
}

void
points_init(points_t p)
{
  p->thetas = new mpfr_t[(p->boundary + p->interior)];
  p->phis = new mpfr_t[(p->boundary + p->interior)];
  for (int i = 0; i < p->interior + p->boundary; i++) {
    mpfr_init(p->thetas[i]);
    mpfr_init(p->phis[i]);
  }
}

void
points_clear(points_t p)
{
  for (int i = 0; i < p->interior + p->boundary; i++) {
    mpfr_clear(p->thetas[i]);
    mpfr_clear(p->phis[i]);
  }
  delete [] p->thetas;
  delete [] p->phis;
}

void
boundary(points_t p, geom_t g)
{
  mpfr_t arcx, arcy, arcz, t, norm, tmp;
  int n;

  mpfr_init(arcx);
  mpfr_init(arcy);
  mpfr_init(arcz);
  mpfr_init(t);
  mpfr_init(norm);
  mpfr_init(tmp);

  n = p->boundary;

  for (int i = 0; i < n; i++) {
    mpfr_set_si(t, i + 1, MPFR_RNDN);
    if (g->half_boundary)
      mpfr_div_si(t, t, 2*n, MPFR_RNDN);
    else
      mpfr_div_si(t, t, n + 1, MPFR_RNDN);

    mpfr_sub(arcx, g->v3[0], g->v2[0], MPFR_RNDN);
    mpfr_mul(arcx, arcx, t, MPFR_RNDN);
    mpfr_add(arcx, arcx, g->v2[0], MPFR_RNDN);

    mpfr_sub(arcy, g->v3[1], g->v2[1], MPFR_RNDN);
    mpfr_mul(arcy, arcy, t, MPFR_RNDN);
    mpfr_add(arcy, arcy, g->v2[1], MPFR_RNDN);

    mpfr_sub(arcz, g->v3[2], g->v2[2], MPFR_RNDN);
    mpfr_mul(arcz, arcz, t, MPFR_RNDN);
    mpfr_add(arcz, arcz, g->v2[2], MPFR_RNDN);

    mpfr_sqr(norm, arcx, MPFR_RNDN);

    mpfr_sqr(tmp, arcy, MPFR_RNDN);
    mpfr_add(norm, norm, tmp, MPFR_RNDN);

    mpfr_sqr(tmp, arcz, MPFR_RNDN);
    mpfr_add(norm, norm, tmp, MPFR_RNDN);

    mpfr_sqrt(norm, norm, MPFR_RNDN);

    mpfr_div(arcx, arcx, norm, MPFR_RNDN);
    mpfr_div(arcy, arcy, norm, MPFR_RNDN);
    mpfr_div(arcz, arcz, norm, MPFR_RNDN);

    mpfr_acos(p->thetas[i], arcz, MPFR_RNDN);
    mpfr_atan2(p->phis[i], arcy, arcx, MPFR_RNDN);
  }

  mpfr_clear(arcx);
  mpfr_clear(arcy);
  mpfr_clear(arcz);
  mpfr_clear(t);
  mpfr_clear(norm);
  mpfr_clear(tmp);
}

void
interior(points_t p, geom_t g)
{
  mpfr_t s, t, x, y, z, norm, tmp;

  mpfr_init(x);
  mpfr_init(y);
  mpfr_init(z);
  mpfr_init(s);
  mpfr_init(t);
  mpfr_init(norm);
  mpfr_init(tmp);

  for (int i = p->boundary; i < p->boundary + p->interior; i++)
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

        mpfr_sub(x, g->v2[0], g->v1[0], MPFR_RNDN);
        mpfr_mul(x, x, s, MPFR_RNDN);
        mpfr_sub(tmp, g->v3[0], g->v1[0], MPFR_RNDN);
        mpfr_mul(tmp, tmp, t, MPFR_RNDN);
        mpfr_add(x, x, tmp, MPFR_RNDN);
        mpfr_add(x, x, g->v1[0], MPFR_RNDN);

        mpfr_sub(y, g->v2[1], g->v1[1], MPFR_RNDN);
        mpfr_mul(y, y, s, MPFR_RNDN);
        mpfr_sub(tmp, g->v3[1], g->v1[1], MPFR_RNDN);
        mpfr_mul(tmp, tmp, t, MPFR_RNDN);
        mpfr_add(y, y, tmp, MPFR_RNDN);
        mpfr_add(y, y, g->v1[1], MPFR_RNDN);

        mpfr_sub(z, g->v2[2], g->v1[2], MPFR_RNDN);
        mpfr_mul(z, z, s, MPFR_RNDN);
        mpfr_sub(tmp, g->v3[2], g->v1[2], MPFR_RNDN);
        mpfr_mul(tmp, tmp, t, MPFR_RNDN);
        mpfr_add(z, z, tmp, MPFR_RNDN);
        mpfr_add(z, z, g->v1[2], MPFR_RNDN);

        mpfr_sqr(norm, x, MPFR_RNDN);

        mpfr_sqr(tmp, y, MPFR_RNDN);
        mpfr_add(norm, norm, tmp, MPFR_RNDN);

        mpfr_sqr(tmp, z, MPFR_RNDN);
        mpfr_add(norm, norm, tmp, MPFR_RNDN);

        mpfr_sqrt(norm, norm, MPFR_RNDN);

        mpfr_div(x, x, norm, MPFR_RNDN);
        mpfr_div(y, y, norm, MPFR_RNDN);
        mpfr_div(z, z, norm, MPFR_RNDN);

        mpfr_acos(p->thetas[i], z, MPFR_RNDN);
        mpfr_atan2(p->phis[i], y, x, MPFR_RNDN);
    }

  mpfr_clear(x);
  mpfr_clear(y);
  mpfr_clear(z);
  mpfr_clear(s);
  mpfr_clear(t);
  mpfr_clear(norm);
  mpfr_clear(tmp);
}
