#ifndef GEOM
#define GEOM

#include "arb.h"

typedef struct
{
  arb_ptr v1;
  arb_ptr v2;
  arb_ptr v3;
  arb_t theta_bound;
  fmpq * angles;
  int half_boundary;
} geom_struct;

typedef geom_struct geom_t[1];

typedef struct
{
  arb_ptr thetas;
  arb_ptr phis;
  int boundary;
  int interior;
} points_struct;

typedef points_struct points_t[1];

void geom_init(geom_t g);

void geom_clear(geom_t g);

void geom_set_angles(geom_t g, int angles_coefs[]);

void geom_compute(geom_t g, slong prec);

void points_init(points_t p, int boundary, int interior);

void points_clear(points_t p);

void boundary(points_t p, geom_t g, slong prec);

void interior(points_t p, geom_t g, slong prec);

#endif
