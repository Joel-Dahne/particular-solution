#ifndef GEOM
#define GEOM

#include "mpfr.h"
#include "stdlib.h"

typedef struct
{
    mpfr_t *v1;
    mpfr_t *v2;
    mpfr_t *v3;
    mpfr_t theta_bound;
    int angles[6];
    int half_boundary;
} geom_struct;

typedef geom_struct geom_t[1];

typedef struct
{
    mpfr_t *thetas;
    mpfr_t *phis;
    int boundary;
    int interior;
} points_struct;

typedef points_struct points_t[1];

void geom_init(geom_t g);

void geom_clear(geom_t g);

void geom_set_prec(geom_t g, mpfr_prec_t prec);

void geom_set(geom_t g, int angles_coefs[]);

void points_init(points_t p, int boundary, int interior);

void points_clear(points_t p);

void boundary(points_t p, geom_t g);

void interior(points_t p, geom_t g);

#endif
