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
 * arb_ptr v1[3]: Coordinates for the vertex placed on the north pole
 * (this will always be [0, 0, 1]).
 *
 * arb_ptr v2[3]: Coordinates for the vertex having y equal to zero.
 *
 * arb_ptr v3[3]: Coordinates for the third vertex.
 *
 * arb_ptr theta_lower[3]: Lower bound of theta for the boundary of
 * the triangle between v2 and v3.
 *
 * arb_ptr theta_upper[3]: Upper bound of theta for the boundary of
 * the triangle between v2 and v3.
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
  arb_ptr theta_lower[3];
  arb_ptr theta_upper[3];
  int vertices[3] = {0};
  int half_edge[3] = {0};
} geom_struct;

typedef geom_struct geom_t[1];

/* The points_t type holds points on the boundary or in the interior
 * of the spherical triangle.
 *
 * In the computations, expansions from one, two or all three vertices
 * are used. Each expansion uses a different parameterization of the
 * sphere (the north pole of the sphere is placed at the vertex from
 * which the expansion originates and one of the other vertices is
 * taken to have phi equal to zero). Because of this all value which
 * depend on the parameterization has three instances, one for each
 * choice of vertex to place at the north pole.
 *
 * Not all points are relevant for all computations. When using an
 * expansion from a particular vertex only the points on the opposite
 * edge and in the interior are relevant, for the other two boundaries
 * the function used will be identically equal to zero. To avoid
 * unnecessary computations we, for each different parameterization,
 * set the value of all points known to evaluate to zero to an
 * indeterminate ball. We then expect the function that uses these
 * values to know that that means that they should evaluate to zero.
 *
 * slong boundary: The total number of points for each edge included
 * in the expansion, that is each edge for which an expansion from the
 * opposite vertex is used. The same number of points are used for
 * each edge and can thus be computed by dividing this number by the
 * number of included edges.
 *
 * slong interior: The total number of interior points.
 *
 * slong total: The total number of interior points.
 *
 * arb_ptr thetas, phi: Holds theta and phi values for all the stored
 * points in the three different parameterizations. Let b be the
 * number of boundary points and i the number of interior points. Then
 * the first b points are on the first edge (the one between the first
 * and second vertex), the b points after that are on the second edge
 * (the one between the second and third vertex) and the b points
 * after that on the third edge (the one between the third and first
 * vertex). After that comes i points from the interior.
*/

typedef struct
{
  slong boundary;
  slong interior;
  slong total;
  arb_ptr thetas[3];
  arb_ptr phis[3];
} points_struct;

typedef points_struct points_t[1];

void geom_init(geom_t g);

void geom_clear(geom_t g);

void geom_set_angles(geom_t g, slong angles_coefs[]);

void geom_compute(geom_t g, slong prec);

void geom_get_mu(arb_t mu, geom_t g, slong vertex, slong i, slong prec);

void points_init(points_t p, geom_t g, slong boundary, slong interior);

void points_clear(points_t p);

void boundary(points_t p, geom_t g, slong prec);

void interior(points_t p, geom_t g, slong prec);

#endif
