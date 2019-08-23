#include "geom.h"

void
geom_init(geom_t g)
{
  g->angle = _fmpq_vec_init(1);
  g->side_length = _arb_vec_init(1);
  g->x = _arb_vec_init(1);
  g->y = _arb_vec_init(1);
  g->r_lower = _arb_vec_init(1);
}

void
geom_clear(geom_t g)
{
  _fmpq_vec_clear(g->angle, 1);
  _arb_vec_clear(g->side_length, 1);
  _arb_vec_clear(g->x, 1);
  _arb_vec_clear(g->y, 1);
  _arb_vec_clear(g->r_lower, 1);
}

/* In practice we only set one angle but keep is like this for
 * compatibility */
void
geom_set_angles(geom_t g, slong angles[])
{
  fmpq_set_si(g->angle, angles[0], angles[1]);
}

void
geom_compute(geom_t g, slong prec)
{
  arb_t t, tmp;
  arb_init(t);
  arb_init(tmp);

  /* Compute x */
  arb_set_fmpq(g->x, g->angle, prec);
  arb_cos_pi(g->x, g->x, prec);
  arb_mul(g->x, g->x, g->side_length, prec);

  /* Compute y */
  arb_set_fmpq(g->y, g->angle, prec);
  arb_sin_pi(g->y, g->y, prec);
  arb_mul(g->y, g->y, g->side_length, prec);

  /* Compute lower bound of r */
  /* With u = t*[x, y] + (1 - t)*[1, 0] compute the critical point
   * with respect to t. */
  /* t = x/((x - 1)^2 + y^2)*/
  arb_sub_si(t, g->x, 1, prec);
  arb_sqr(t, t, prec);
  arb_addmul(t, g->y, g->y, prec);
  arb_div(t, g->x, t, prec);

  /* If the critical point is on the interval (0, 1) the lower bound
   * of r is given by the value at this point. Otherwise it's given by
   * the smaller of the two side lengths. */
  arb_sub_si(tmp, t, 1, prec);
  if (arb_contains_positive(t) && arb_contains_negative(tmp))
  {
    arb_mul(g->r_lower, t, g->x, prec);
    arb_add_si(g->r_lower, g->r_lower, 1, prec);
    arb_sub(g->r_lower, g->r_lower, t, prec);
    arb_sqr(g->r_lower, g->r_lower, prec);
    arb_mul(tmp, t, g->y, prec);
    arb_sqr(tmp, tmp, prec);
    arb_add(g->r_lower, g->r_lower, tmp, prec);
  }
  else{
    arb_set_si(tmp, 1);
    arb_min(g->r_lower, g->side_length, tmp, prec);
  }

  arb_clear(t);
  arb_clear(tmp);
}

/* mu would in the plane normally be called nu */
void
geom_get_mu(arb_t nu, geom_t g, slong vertex, slong i, slong prec)
{
  fmpq_t tmp;

  fmpq_init(tmp);

  fmpq_set_si(tmp, i + 1, 1);
  fmpq_div(tmp, tmp, g->angle);
  arb_set_fmpq(nu, tmp, prec);

  fmpq_clear(tmp);
}

void
geom_area(arb_t area, geom_t g, slong prec)
{
  arb_div_si(area, g->y, 2, prec);
}

void
points_init(points_t p, geom_t g, slong boundary, slong interior)
{
  p->boundary = boundary;
  p->interior = interior;
  p->total = p->boundary + p->interior;
  for (slong i = 0; i < 3; i++)
  {
    p->thetas[i] = _arb_vec_init(p->total);
    p->phis[i] = _arb_vec_init(p->total);
  }
}

void
points_clear(points_t p)
{
  for (slong i = 0; i < 3; i++)
  {
    _arb_vec_clear(p->thetas[i], p->total);
    _arb_vec_clear(p->phis[i], p->total);
  }
}

/* Boundary and interior points are at the moment only implemented in
 * Julia. */
void
boundary(points_t p, geom_t g, slong prec)
{

}

/* Boundary and interior points are at the moment only implemented in
 * Julia. */
void
interior(points_t p, geom_t g, slong prec)
{

}

/* TODO: Implement this */
void
parametrization(arb_ptr r, arb_ptr theta, geom_t geom, arb_t t, slong n,
                slong vertex, slong prec)
{

}
