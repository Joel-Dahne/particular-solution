#ifndef GEOM
#define GEOM

#include "arb.h"

/* The triangle is given by putting one vertex at the origin, one at
 * (1, 0) and the third vertex is placed so that the angle at the
 * origin is equal to the given one and the distance between the
 * origin and the new vertex is the given side_length. The variables x
 * and y are computed to correspond to this vertex. For compatibility
 * with the spherical triangles we need to keep the vertices and
 * half_edge variables but don't use them in practice. */
typedef struct
{
  fmpq * angle;
  arb_ptr side_length;
  arb_ptr x;
  arb_ptr y;
  int vertices[3];
  int half_edge[3];
} geom_struct;

typedef geom_struct geom_t[1];

typedef struct
{
  slong boundary;
  slong interior;
  slong total;
  arb_ptr thetas[3]; // radiis
  arb_ptr phis[3]; // angles
} points_struct;

typedef points_struct points_t[1];

void geom_init(geom_t g);

void geom_clear(geom_t g);

void geom_set_angles(geom_t g, slong angles_coefs[]);

void geom_compute(geom_t g, slong prec);

void geom_get_mu(arb_t nu, geom_t g, slong vertex, slong i, slong prec);

void geom_area(arb_t area, geom_t g, slong prec);

void points_init(points_t p, geom_t g, slong boundary, slong interior);

void points_clear(points_t p);

void boundary(points_t p, geom_t g, slong prec);

void interior(points_t p, geom_t g, slong prec);

void parametrization(arb_ptr r, arb_ptr theta, geom_t geom, arb_t t, slong n,
                     slong vertex, slong prec);

#endif
