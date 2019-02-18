#include "geom.h"

void
geom_init(geom_t g)
{
  g->angles = _fmpq_vec_init(3);
  for (slong i = 0; i < 3; i++)
  {
    g->v1[i] = _arb_vec_init(3);
    g->v2[i] = _arb_vec_init(3);
    g->v3[i] = _arb_vec_init(3);
  }
}

void
geom_clear(geom_t g)
{
  _fmpq_vec_clear(g->angles, 3);
    for (slong i = 0; i < 3; i++)
  {
    _arb_vec_clear(g->v1[i], 3);
    _arb_vec_clear(g->v2[i], 3);
    _arb_vec_clear(g->v3[i], 3);
  }
}

void
geom_set_angles(geom_t g, slong angles[])
{
  fmpq_set_si(g->angles + 0, angles[0], angles[1]);
  fmpq_set_si(g->angles + 1, angles[2], angles[3]);
  fmpq_set_si(g->angles + 2, angles[4], angles[5]);
}

void
geom_compute(geom_t g, slong prec)
{
  arb_ptr angles;
  arb_t S, tmp1, tmp2;

  angles = _arb_vec_init(3);

  arb_init(S);
  arb_init(tmp1);
  arb_init(tmp2);

  for (slong i = 0; i < 3; i++)
  {
    /* Compute the angles from the quotient */
    arb_const_pi(angles + 0, prec);
    arb_const_pi(angles + 1, prec);
    arb_const_pi(angles + 2, prec);

    /* Angle for the vertex at the north pole in the current
     * parameterization. */
    arb_mul_fmpz(angles + 0, angles + 0, fmpq_numref(g->angles + (3 - i)%3), prec);
    arb_div_fmpz(angles + 0, angles + 0, fmpq_denref(g->angles + (3 - i)%3), prec);
    /* Angle for the vertex having y equal to zero in the current
     * parameterization */
    arb_mul_fmpz(angles + 1, angles + 1, fmpq_numref(g->angles + (4 - i)%3), prec);
    arb_div_fmpz(angles + 1, angles + 1, fmpq_denref(g->angles + (4 - i)%3), prec);
    /* Angle for the third vertex in the current parameterization */
    arb_mul_fmpz(angles + 2, angles + 2, fmpq_numref(g->angles + (5 - i)%3), prec);
    arb_div_fmpz(angles + 2, angles + 2, fmpq_denref(g->angles + (5 - i)%3), prec);

    /* Compute the vectors for the spherical triangle */
    /* S = (angles[0] + angles[1] + angles[2])/2 */
    arb_add(S, angles + 0, angles + 1, prec);
    arb_add(S, S, angles + 2, prec);
    arb_div_si(S, S, 2, prec);

    /* v1 = [0, 0, 1] */
    arb_set_si(g->v1[i] + 0, 0);
    arb_set_si(g->v1[i] + 1, 0);
    arb_set_si(g->v1[i] + 2, 1);

    /* Compute theta for v2 */
    arb_cos(tmp1, S, prec);

    arb_sub(tmp2, S, angles + 1, prec);
    arb_cos(tmp2, tmp2, prec);
    arb_mul(tmp1, tmp1, tmp2, prec);

    arb_sin(tmp2, angles + 0, prec);
    arb_div(tmp1, tmp1, tmp2, prec);

    arb_sin(tmp2, angles + 2, prec);
    arb_div(tmp1, tmp1, tmp2, prec);

    arb_neg(tmp1, tmp1);
    arb_sqrt(tmp1, tmp1, prec);
    arb_asin(tmp1, tmp1, prec);
    arb_mul_si(tmp1, tmp1, 2, prec);

    /* Compute v2 from knowing theta */
    arb_sin(g->v2[i] + 0, tmp1, prec);
    arb_set_si(g->v2[i] + 1, 0);
    arb_cos(g->v2[i] + 2, tmp1, prec);

    /* Compute theta for v3 */
    arb_cos(tmp1, S, prec);

    arb_sub(tmp2, S, angles + 2, prec);
    arb_cos(tmp2, tmp2, prec);
    arb_mul(tmp1, tmp1, tmp2, prec);

    arb_sin(tmp2, angles + 0, prec);
    arb_div(tmp1, tmp1, tmp2, prec);

    arb_sin(tmp2, angles + 1, prec);
    arb_div(tmp1, tmp1, tmp2, prec);

    arb_neg(tmp1, tmp1);
    arb_sqrt(tmp1, tmp1, prec);
    arb_asin(tmp1, tmp1, prec);
    arb_mul_si(tmp1, tmp1, 2, prec);

    /* Compute v3 from knowing theta */
    arb_sin(tmp2, tmp1, prec);
    arb_cos(g->v3[i] + 0, angles + 0, prec);
    arb_mul(g->v3[i] + 0, g->v3[i] + 0, tmp2, prec);
    arb_sin(g->v3[i] + 1, angles + 0, prec);
    arb_mul(g->v3[i] + 1, g->v3[i] + 1, tmp2, prec);
    arb_cos(g->v3[i] + 2, tmp1, prec);
  }

  arb_clear(S);
  arb_clear(tmp1);
  arb_clear(tmp2);

  _arb_vec_clear(angles, 3);
}

void
geom_get_mu(arb_t mu, geom_t g, slong vertex, slong i, slong prec)
{
  fmpq_t tmp;

  fmpq_init(tmp);

  fmpq_set_si(tmp, -(g->half_edge[vertex] ? 2*i + 1 : i + 1), 1);
  fmpq_div(tmp, tmp, g->angles + (3 - vertex)%3);
  arb_set_fmpq(mu, tmp, prec);

  fmpq_clear(tmp);
}

void
points_init(points_t p, slong boundary, slong interior)
{
  p->boundary = boundary;
  p->interior = interior;
  p->thetas = _arb_vec_init(p->boundary + p->interior);
  p->phis = _arb_vec_init(p->boundary + p->interior);
}

void
points_clear(points_t p)
{
  _arb_vec_clear(p->thetas, p->interior + p->boundary);
  _arb_vec_clear(p->phis, p->interior + p->boundary);
}

void
boundary(points_t p, geom_t g, slong prec)
{
  arb_ptr arc;
  arb_t t, norm, tmp;
  int n;

  arc = _arb_vec_init(3);

  arb_init(t);
  arb_init(norm);
  arb_init(tmp);

  n = p->boundary;

  for (int i = 0; i < n; i++) {
    arb_set_si(t, i + 1);
    if (g->half_edge[0])
      arb_div_si(t, t, 2*n, prec);
    else
      arb_div_si(t, t, n + 1, prec);

    _arb_vec_sub(arc, g->v3[0], g->v2[0], 3, prec);
    _arb_vec_scalar_mul(arc, arc, 3, t, prec);
    _arb_vec_add(arc, arc, g->v2[0], 3, prec);

    _arb_vec_dot(norm, arc, arc, 3, prec);
    arb_sqrt(norm, norm, prec);

    _arb_vec_scalar_div(arc, arc, 3, norm, prec);

    arb_acos(p->thetas + i, arc + 2, prec);
    arb_atan2(p->phis + i, arc + 1, arc + 0, prec);
  }

  _arb_vec_clear(arc, 3);

  arb_clear(t);
  arb_clear(norm);
  arb_clear(tmp);
}

void
interior(points_t p, geom_t g, slong prec)
{
  arb_ptr xyz, v2mv1, v3mv1;
  arb_t s, t, norm, tmp;

  xyz = _arb_vec_init(3);
  v2mv1 = _arb_vec_init(3);
  v3mv1 = _arb_vec_init(3);

  arb_init(s);
  arb_init(t);
  arb_init(norm);
  arb_init(tmp);

  _arb_vec_sub(v2mv1, g->v2[0], g->v1[0], 3, prec);
  _arb_vec_sub(v3mv1, g->v3[0], g->v1[0], 3, prec);

  for (int i = p->boundary; i < p->boundary + p->interior; i++)
  {
    /* We take s and t random with s in [0, 1) and s < 1 - t*/
    arb_set_d(s, ((double)rand()/RAND_MAX));
    arb_set_d(t, ((double)rand()/RAND_MAX));

    /* Set tmp to 1 - t */
    arb_sub_si(tmp, t, 1, prec);
    arb_neg(tmp, tmp);
    /* If s >= 1 - t set s = 1 - t, t = 1 - s to make them satisfy
       s < 1 - t. */
    if (arb_ge(s, tmp))
    {
      /* Set t to 1 - s */
      arb_sub_si(t, s, 1, prec);
      arb_neg(t, t);

      /* Set s to 1 - t */
      arb_set(s, tmp);
    }

    _arb_vec_set(xyz, g->v1[0], 3);
    _arb_vec_scalar_addmul(xyz, v2mv1, 3, s, prec);
    _arb_vec_scalar_addmul(xyz, v3mv1, 3, t, prec);

    _arb_vec_dot(norm, xyz, xyz, 3, prec);
    arb_sqrt(norm, norm, prec);

    _arb_vec_scalar_div(xyz, xyz, 3, norm, prec);

    arb_acos(p->thetas + i, xyz + 2, prec);
    arb_atan2(p->phis + i, xyz + 1, xyz + 0, prec);
  }

  _arb_vec_clear(xyz, 3);
  _arb_vec_clear(v2mv1, 3);
  _arb_vec_clear(v3mv1, 3);

  arb_clear(s);
  arb_clear(t);
  arb_clear(norm);
  arb_clear(tmp);
}
