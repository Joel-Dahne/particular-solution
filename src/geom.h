#ifndef GEOM
#define GEOM

#include "arb.h"

/* The geom_t type represents the geometry of the spherical triangle.
 * The triangle is defined by the angles at the three vertices.
 *
 * In the computations, expansions from one, two or all three vertices
 * are used. Each expansion uses a different parameterization of the
 * sphere (the north pole of the sphere is placed at the vertex from
 * which the expansion originates and one of the other vertices is
 * taken to have phi equal to zero). Because of this all value which
 * depend on the parameterization has three instances, one for each
 * choice of vertex to place at the north pole.
 *
 * fmpq * angles: Contains the three angles for the three different
 * vertices.
 *
 * arb_ptr v1[3]: Coordinates for the vertex places on the north pole
 * (this will always be [0, 0, 1]).
 *
 * arb_ptr v2[3]: Coordinates for the vertex having y equal to zero.
 *
 * arb_ptr v3[3]: Coordinates for the third vertex.
 *
 * int vertices[3]: Flags for which vertices to use expansions from (1
 * to use it, 0 to not use it).
 *
 * int half_edge[3]: Flags for which expansions only half the the edge
 * is to be used in the computations, in general because of symmetry
 * reasons.
 */

typedef struct
{
  fmpq * angles;
  arb_ptr v1[3];
  arb_ptr v2[3];
  arb_ptr v3[3];
  int vertices[3] = {0};
  int half_edge[3] = {0};
} geom_struct;

typedef geom_struct geom_t[1];

typedef struct
{
  arb_ptr thetas;
  arb_ptr phis;
  slong boundary;
  slong interior;
} points_struct;

typedef points_struct points_t[1];

void geom_init(geom_t g);

void geom_clear(geom_t g);

void geom_set_angles(geom_t g, slong angles_coefs[]);

void geom_compute(geom_t g, slong prec);

void points_init(points_t p, slong boundary, slong interior);

void points_clear(points_t p);

void boundary(points_t p, geom_t g, slong prec);

void interior(points_t p, geom_t g, slong prec);

#endif
