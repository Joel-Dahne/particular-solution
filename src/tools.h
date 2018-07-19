#ifndef TOOLS
#define TOOLS

#include "mpfr.h"

struct Geometry {
    mpfr_t v1[3];
    mpfr_t v2[3];
    mpfr_t v3[3];
    mpfr_t theta_bound;
    int half_boundary;
};

struct Points {
    mpfr_t *thetas;
    mpfr_t *phis;
    int boundary;
    int interior;
};

void points_init(struct Points & points);

void points_clear(struct Points & points);

void geometry_init(struct Geometry & geometry);

void geometry_set_prec(struct Geometry & geometry, mpfr_prec_t prec);

void geometry_clear(struct Geometry & geometry);

void angles_to_vectors(struct Geometry &geometry, mpfr_t *angles);

void boundary(struct Points points, struct Geometry geometry);

void interior(struct Points points, struct Geometry geometry);

#endif
