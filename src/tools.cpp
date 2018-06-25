#include "mpfr.h"
#include "stdlib.h"

/* The vector angles should contain three angles and the vectors v1,
 * v2, v3 are set to the coordinates of the vertices of the spherical
 * triangle having the angles defined in angles. One vertex is chosen
 * to be the north pole, one to have phi=0 and the third to have
 * phi=angles[0]. The parameter theta_bound is set to the average of
 * the theta values for v2 and v3, it is mainly intended to be used
 * for computing the scaling integrals. The vectors v1, v2, v3 should
 * be initialized to have room for 3 elements. */
void
angles_to_vectors(mpfr_t *v1, mpfr_t *v2, mpfr_t *v3, mpfr_t theta_bound,
                  mpfr_t *angles)
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
  mpfr_set_si(v1[0], 0, MPFR_RNDN);
  mpfr_set_si(v1[1], 0, MPFR_RNDN);
  mpfr_set_si(v1[2], 1, MPFR_RNDN);

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
  mpfr_set(theta_bound, tmp1, MPFR_RNDN);

  /* Compute v2 from knowing theta */
  mpfr_sin(v2[0], tmp1, MPFR_RNDN);
  mpfr_set_si(v2[1], 0, MPFR_RNDN);
  mpfr_cos(v2[2], tmp1, MPFR_RNDN);

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
  mpfr_add(theta_bound, theta_bound, tmp1, MPFR_RNDN);
  mpfr_div_si(theta_bound, theta_bound, 2, MPFR_RNDN);

  /* Compute v3 from knowing theta */
  mpfr_sin(tmp2, tmp1, MPFR_RNDN);
  mpfr_cos(v3[0], angles[0], MPFR_RNDN);
  mpfr_mul(v3[0], v3[0], tmp2, MPFR_RNDN);
  mpfr_sin(v3[1], angles[0], MPFR_RNDN);
  mpfr_mul(v3[1], v3[1], tmp2, MPFR_RNDN);
  mpfr_cos(v3[2], tmp1, MPFR_RNDN);

  mpfr_clear(S);
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);
}

void boundary(mpfr_t *thetas, mpfr_t *phis, mpfr_t *v1, mpfr_t *v2, int n,
              int half_boundary) {
  mpfr_t arcx, arcy, arcz, t, norm, tmp;

  mpfr_init(arcx);
  mpfr_init(arcy);
  mpfr_init(arcz);

  mpfr_init(t);
  mpfr_init(norm);
  mpfr_init(tmp);

  for (int i = 0; i < n; i++) {
    mpfr_set_si(t, i + 1, MPFR_RNDN);
    if (half_boundary)
      mpfr_div_si(t, t, 2*n, MPFR_RNDN);
    else
      mpfr_div_si(t, t, n + 1, MPFR_RNDN);

    mpfr_sub(arcx, v2[0], v1[0], MPFR_RNDN);
    mpfr_mul(arcx, arcx, t, MPFR_RNDN);
    mpfr_add(arcx, arcx, v1[0], MPFR_RNDN);

    mpfr_sub(arcy, v2[1], v1[1], MPFR_RNDN);
    mpfr_mul(arcy, arcy, t, MPFR_RNDN);
    mpfr_add(arcy, arcy, v1[1], MPFR_RNDN);

    mpfr_sub(arcz, v2[2], v1[2], MPFR_RNDN);
    mpfr_mul(arcz, arcz, t, MPFR_RNDN);
    mpfr_add(arcz, arcz, v1[2], MPFR_RNDN);

    mpfr_sqr(norm, arcx, MPFR_RNDN);

    mpfr_sqr(tmp, arcy, MPFR_RNDN);
    mpfr_add(norm, norm, tmp, MPFR_RNDN);

    mpfr_sqr(tmp, arcz, MPFR_RNDN);
    mpfr_add(norm, norm, tmp, MPFR_RNDN);

    mpfr_sqrt(norm, norm, MPFR_RNDN);

    mpfr_div(arcx, arcx, norm, MPFR_RNDN);
    mpfr_div(arcy, arcy, norm, MPFR_RNDN);
    mpfr_div(arcz, arcz, norm, MPFR_RNDN);

    mpfr_acos(thetas[i], arcz, MPFR_RNDN);
    mpfr_atan2(phis[i], arcy, arcx, MPFR_RNDN);
  }

  mpfr_clear(arcx);
  mpfr_clear(arcy);
  mpfr_clear(arcz);

  mpfr_clear(t);
  mpfr_clear(norm);
  mpfr_clear(tmp);
}

void interior(mpfr_t *thetas, mpfr_t *phis, mpfr_t *v1, mpfr_t *v2, mpfr_t *v3,
              int n) {
  mpfr_t s, t, x, y, z, norm, tmp;

  mpfr_init(x);
  mpfr_init(y);
  mpfr_init(z);

  mpfr_init(s);
  mpfr_init(t);
  mpfr_init(norm);
  mpfr_init(tmp);

  for (int i = 0; i < n; i++)
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

        mpfr_sub(x, v2[0], v1[0], MPFR_RNDN);
        mpfr_mul(x, x, s, MPFR_RNDN);
        mpfr_sub(tmp, v3[0], v1[0], MPFR_RNDN);
        mpfr_mul(tmp, tmp, t, MPFR_RNDN);
        mpfr_add(x, x, tmp, MPFR_RNDN);
        mpfr_add(x, x, v1[0], MPFR_RNDN);

        mpfr_sub(y, v2[1], v1[1], MPFR_RNDN);
        mpfr_mul(y, y, s, MPFR_RNDN);
        mpfr_sub(tmp, v3[1], v1[1], MPFR_RNDN);
        mpfr_mul(tmp, tmp, t, MPFR_RNDN);
        mpfr_add(y, y, tmp, MPFR_RNDN);
        mpfr_add(y, y, v1[1], MPFR_RNDN);

        mpfr_sub(z, v2[2], v1[2], MPFR_RNDN);
        mpfr_mul(z, z, s, MPFR_RNDN);
        mpfr_sub(tmp, v3[2], v1[2], MPFR_RNDN);
        mpfr_mul(tmp, tmp, t, MPFR_RNDN);
        mpfr_add(z, z, tmp, MPFR_RNDN);
        mpfr_add(z, z, v1[2], MPFR_RNDN);

        mpfr_sqr(norm, x, MPFR_RNDN);

        mpfr_sqr(tmp, y, MPFR_RNDN);
        mpfr_add(norm, norm, tmp, MPFR_RNDN);

        mpfr_sqr(tmp, z, MPFR_RNDN);
        mpfr_add(norm, norm, tmp, MPFR_RNDN);

        mpfr_sqrt(norm, norm, MPFR_RNDN);

        mpfr_div(x, x, norm, MPFR_RNDN);
        mpfr_div(y, y, norm, MPFR_RNDN);
        mpfr_div(z, z, norm, MPFR_RNDN);

        mpfr_acos(thetas[i], z, MPFR_RNDN);
        mpfr_atan2(phis[i], y, x, MPFR_RNDN);
    }

  mpfr_clear(x);
  mpfr_clear(y);
  mpfr_clear(z);

  mpfr_clear(s);
  mpfr_clear(t);
  mpfr_clear(norm);
  mpfr_clear(tmp);
}
