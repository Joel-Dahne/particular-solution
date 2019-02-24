#include "geom.h"
#include "arb_mat.h"

void
geom_init(geom_t g)
{
  g->angles = _fmpq_vec_init(3);
  for (slong i = 0; i < 3; i++)
  {
    g->v1[i] = _arb_vec_init(3);
    g->v2[i] = _arb_vec_init(3);
    g->v3[i] = _arb_vec_init(3);
    g->theta_lower[i] = _arb_vec_init(3);
    g->theta_upper[i] = _arb_vec_init(3);
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
    _arb_vec_clear(g->theta_lower[i], 3);
    _arb_vec_clear(g->theta_upper[i], 3);
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
  arb_ptr angles, w;
  arb_t S, theta_v2, theta_v3, critical_point, tmp, tmp_dot;

  angles = _arb_vec_init(3);
  w = _arb_vec_init(3);

  arb_init(S);
  arb_init(theta_v2);
  arb_init(theta_v3);
  arb_init(critical_point);
  arb_init(tmp);
  arb_init(tmp_dot);

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

    /* Compute the theta value for for v2 */
    arb_cos(theta_v2, S, prec);

    arb_sub(tmp, S, angles + 1, prec);
    arb_cos(tmp, tmp, prec);
    arb_mul(theta_v2, theta_v2, tmp, prec);

    arb_sin(tmp, angles + 0, prec);
    arb_div(theta_v2, theta_v2, tmp, prec);

    arb_sin(tmp, angles + 2, prec);
    arb_div(theta_v2, theta_v2, tmp, prec);

    arb_neg(theta_v2, theta_v2);
    arb_sqrt(theta_v2, theta_v2, prec);
    arb_asin(theta_v2, theta_v2, prec);
    arb_mul_si(theta_v2, theta_v2, 2, prec);

    /* Compute v2 from knowing theta
    * v2 = [sin(theta_v2), 0, cos(theta_v2)] */
    arb_sin(g->v2[i] + 0, theta_v2, prec);
    arb_set_si(g->v2[i] + 1, 0);
    arb_cos(g->v2[i] + 2, theta_v2, prec);

    /* Compute the theta value for v3 */
    arb_cos(theta_v3, S, prec);

    arb_sub(tmp, S, angles + 2, prec);
    arb_cos(tmp, tmp, prec);
    arb_mul(theta_v3, theta_v3, tmp, prec);

    arb_sin(tmp, angles + 0, prec);
    arb_div(theta_v3, theta_v3, tmp, prec);

    arb_sin(tmp, angles + 1, prec);
    arb_div(theta_v3, theta_v3, tmp, prec);

    arb_neg(theta_v3, theta_v3);
    arb_sqrt(theta_v3, theta_v3, prec);
    arb_asin(theta_v3, theta_v3, prec);
    arb_mul_si(theta_v3, theta_v3, 2, prec);

    /* Compute v3 from knowing theta
    * v3 = [sin(theta_v3)*cos(angles[0]), sin(theta_v3)*sin(angles[0]),
    *       cos(theta_v3)] */
    arb_sin(tmp, theta_v3, prec);
    arb_cos(g->v3[i] + 0, angles + 0, prec);
    arb_mul(g->v3[i] + 0, g->v3[i] + 0, tmp, prec);
    arb_sin(g->v3[i] + 1, angles + 0, prec);
    arb_mul(g->v3[i] + 1, g->v3[i] + 1, tmp, prec);
    arb_cos(g->v3[i] + 2, theta_v3, prec);

    /* Compute the critical point, with respect to theta, on the
     * boundary between v2 and v3. */
    /* w = v3 - v2 */
    _arb_vec_sub(w, g->v3[i], g->v2[i], 3, prec);
    /* dot(v2, v2) */
    _arb_vec_dot(tmp_dot, g->v2[i], g->v2[i], 3, prec);
    /* critical_point = w_3*dot(v2, v2) */
    arb_mul(critical_point, w + 2, tmp_dot, prec);
    /* dot(v2, w) */
    _arb_vec_dot(tmp_dot, g->v2[i], w, 3, prec);
    /* v2_3*dot(v2, w) */
    arb_mul(tmp, g->v2[i] + 2, tmp_dot, prec);
    /* critical_point = v_3*dit(v2, w) - w_3*dot(v2, v2) */
    arb_sub(critical_point, tmp, critical_point, prec);
    /* w_3*dot(v2, w) */
    arb_mul(tmp, w + 2, tmp_dot, prec);
    /* dot(w, w) */
    _arb_vec_dot(tmp_dot, w, w, 3, prec);
    /* w_3*dot(v2, w) - v_3*dot(w, w) */
    arb_submul(tmp, g->v2[i] + 2, tmp_dot, prec);

    /* TODO: Handle the case when the enclosure of tmp contains 0. If
     * tmp is exactly zero we can safetly discard the critical point.
     * However since we only have an enclosure we can't be sure tmp is
     * exactly zero. Most likely it will in theory be exactly zero and
     * there are no solutions, but we need to prove this. */
    /* The critical point */
    arb_div(critical_point, critical_point, tmp, prec);

    /* If the critical point is in the interval (0, 1) compute the angle
     * corresponding to the critical point */
    arb_sub_si(tmp, critical_point, 1, prec);
    if (arb_contains_positive(critical_point) && arb_contains_negative(tmp))
    {
      /* Compute the vector v2 + w*critical_point = v2 + (v3 -
       * v2)*critical point, i.e. the vector at the critical point. */
      _arb_vec_scalar_mul(w, w, 3, critical_point, prec);
      _arb_vec_add(w, g->v2[i], w, 3, prec);
      /* Comute the norm of the above vector */
      _arb_vec_dot(tmp_dot, w, w, 3, prec);
      arb_sqrt(tmp_dot, tmp_dot, prec);
      /* Compute z for this point */
      arb_div(tmp, w + 2, tmp_dot, prec);

      /* Compute the theta value for this point */
      arb_acos(tmp, tmp, prec);
    }
    else
    {
      /* There is no critical point on the edge so put the theta value
       * to that of theta_v2, that way it will not affect the minimum
       * of maximum theta value. */
      arb_set(tmp, theta_v2);
    }

    /* Compute the minimum theta value */
    arb_min(g->theta_lower[i], theta_v2, theta_v3, prec);
    arb_min(g->theta_lower[i], g->theta_lower[i], tmp, prec);
    /* Compute the maximum theta value */
    arb_max(g->theta_upper[i], theta_v2, theta_v3, prec);
    arb_max(g->theta_upper[i], g->theta_upper[i], tmp, prec);
  }

  arb_clear(S);
  arb_clear(theta_v2);
  arb_clear(theta_v3);
  arb_clear(critical_point);
  arb_clear(tmp);
  arb_clear(tmp_dot);

  _arb_vec_clear(angles, 3);
  _arb_vec_clear(w, 3);
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
points_init(points_t p, geom_t g, slong per_boundary, slong interior)
{
  slong boundaries;
  boundaries = g->vertices[0] + g->vertices[1] + g->vertices[2];
  p->boundary = boundaries*per_boundary;
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

/* Computes spherical coordinates (theta, phi) of a point in Cartesian
 * coordinates (x, y, z). First projects the point onto the sphere and
 * then computes the spherical coordinates of that point. */
static void
cartesian_to_spherical(arb_t theta, arb_t phi, arb_ptr xyz, slong prec)
{
  arb_ptr tmp;
  arb_t norm;

  tmp = _arb_vec_init(3);

  arb_init(norm);

  _arb_vec_dot(norm, xyz, xyz, 3, prec);
  arb_sqrt(norm, norm, prec);
  _arb_vec_scalar_div(tmp, xyz, 3, norm, prec);

  arb_acos(theta, tmp + 2, prec);
  arb_atan2(phi, tmp + 1, tmp + 0, prec);

  _arb_vec_clear(tmp, 3);

  arb_clear(norm);
}

void
boundary(points_t p, geom_t g, slong prec)
{
  arb_ptr xyz;
  arb_t t;
  slong boundaries, start;
  int n;

  xyz = _arb_vec_init(3);

  arb_init(t);

  boundaries = g->vertices[0] + g->vertices[1] + g->vertices[2];
  n = p->boundary/boundaries;

  for (slong vertex = 0; vertex < 3; vertex++)
  {
    /* Begin by setting all points to indeterminate balls and then
     * proceed to compute values to the points which are not
     * supposed to be indeterminate. */
    _arb_vec_indeterminate(p->thetas[vertex], p->boundary);
    _arb_vec_indeterminate(p->phis[vertex], p->boundary);

    if (g->vertices[vertex])
    {
      /* Compute starting index for the points which are not supposed to
       * be zero. */
      start = 0;
      for (slong i = 0; i < vertex; i++)
      {
        start += g->vertices[i] ? n : 0;
      }

      for (slong i = start; i < start + n; i++) {
        /* Set value of t */
        arb_set_si(t, i - start + 1);
        if (g->half_edge[vertex])
          arb_div_si(t, t, 2*n, prec);
        else
          arb_div_si(t, t, n + 1, prec);

        /* xyz = v2 - t*(v3 -v2) */
        _arb_vec_sub(xyz, g->v3[vertex], g->v2[vertex], 3, prec);
        _arb_vec_scalar_mul(xyz, xyz, 3, t, prec);
        _arb_vec_add(xyz, xyz, g->v2[vertex], 3, prec);

        cartesian_to_spherical(p->thetas[vertex] + i, p->phis[vertex] + i, xyz,
                               prec);
      }
    }
  }

  _arb_vec_clear(xyz, 3);

  arb_clear(t);
}

void
interior(points_t p, geom_t g, slong prec)
{
  arb_mat_t M[3], xyz_matrix;
  arb_ptr xyz, v2mv1, v3mv1;
  arb_t alpha, beta, cosalpha, cosbeta, sinalpha, sinbeta, s, t, tmp;

  arb_mat_init(M[0], 3, 3);
  arb_mat_init(M[1], 3, 3);
  arb_mat_init(M[2], 3, 3);
  arb_mat_init(xyz_matrix, 3, 1);

  xyz = _arb_vec_init(3);
  v2mv1 = _arb_vec_init(3);
  v3mv1 = _arb_vec_init(3);

  arb_init(alpha);
  arb_init(beta);
  arb_init(cosalpha);
  arb_init(cosbeta);
  arb_init(sinalpha);
  arb_init(sinbeta);
  arb_init(s);
  arb_init(t);
  arb_init(tmp);

  /* We need to compute a number of random points in the interior of
   * the spherical triangle. For each point we need its values in the
   * three different parameterizations.
   *
   * We begin by generating a random point in the Euclidean triangle
   * having vertices in (0, 0), (0, 1) and (1, 0). These points we map
   * to Cartesian coordinates in the first parameterization. We can
   * then map them to Cartesian coordinates in the different
   * parameterizations. Finally we can map these Cartesian coordinates
   * to the required spherical coordinates. */

  _arb_vec_sub(v2mv1, g->v2[0], g->v1[0], 3, prec);
  _arb_vec_sub(v3mv1, g->v3[0], g->v1[0], 3, prec);

  /* Compute matrices for the coordinate change. */

  /* Coordinates 1 to coordinates 1 */
  arb_mat_one(M[0]);

  /* Coordinates 1 to coordinates 2 */
  /* We rotate by alpha along the y-axis and then by beta along the z
   * axis */
  _arb_vec_dot(alpha, g->v1[0], g->v2[0], 3, prec);
  arb_acos(alpha, alpha, prec);
  arb_neg(alpha, alpha);

  arb_set_fmpq(beta, g->angles + 2, prec);
  arb_const_pi(tmp, prec);
  arb_mul(beta, beta, tmp, prec);
  arb_sub(beta, beta, tmp, prec);

  arb_cos(cosalpha, alpha, prec);
  arb_cos(cosbeta, beta, prec);
  arb_sin(sinalpha, alpha, prec);
  arb_sin(sinbeta, beta, prec);

  /* M[1] = [cos(alpha)cos(beta) -sin(beta) cos(beta)sin(alpha)]
   *        [cos(alpha)sin(beta) cos(beta) sin(alpha)sin(beta) ]
   *        [-sin(alpha)          0          cos(alpha)        ] */
  arb_mul(arb_mat_entry(M[1], 0, 0), cosalpha, cosbeta, prec);
  arb_neg(arb_mat_entry(M[1], 0, 1), sinbeta);
  arb_mul(arb_mat_entry(M[1], 0, 2), cosbeta, sinalpha, prec);

  arb_mul(arb_mat_entry(M[1], 1, 0), cosalpha, sinbeta, prec);
  arb_set(arb_mat_entry(M[1], 1, 1), cosbeta);
  arb_mul(arb_mat_entry(M[1], 1, 2), sinalpha, sinbeta, prec);

  arb_neg(arb_mat_entry(M[1], 2, 0), sinalpha);
  arb_zero(arb_mat_entry(M[1], 2, 1));
  arb_set(arb_mat_entry(M[1], 2, 2), cosalpha);

  /* Coordinates 2 to coordinates 3 */
  /* We rotate by alpha along the y-axis and then by beta along the z
   * axis */
  _arb_vec_dot(alpha, g->v1[1], g->v2[1], 3, prec);
  arb_acos(alpha, alpha, prec);
  arb_neg(alpha, alpha);

  arb_set_fmpq(beta, g->angles + 1, prec);
  arb_const_pi(tmp, prec);
  arb_mul(beta, beta, tmp, prec);
  arb_sub(beta, beta, tmp, prec);

  arb_cos(cosalpha, alpha, prec);
  arb_cos(cosbeta, beta, prec);
  arb_sin(sinalpha, alpha, prec);
  arb_sin(sinbeta, beta, prec);

  /* M[2] = [cos(alpha)cos(beta) -sin(beta) cos(beta)sin(alpha)]
   *        [cos(alpha)sin(beta) cos(beta) sin(alpha)sin(beta) ]
   *        [-sin(alpha)          0          cos(alpha)        ] */
  arb_mul(arb_mat_entry(M[2], 0, 0), cosalpha, cosbeta, prec);
  arb_neg(arb_mat_entry(M[2], 0, 1), sinbeta);
  arb_mul(arb_mat_entry(M[2], 0, 2), cosbeta, sinalpha, prec);

  arb_mul(arb_mat_entry(M[2], 1, 0), cosalpha, sinbeta, prec);
  arb_set(arb_mat_entry(M[2], 1, 1), cosbeta);
  arb_mul(arb_mat_entry(M[2], 1, 2), sinalpha, sinbeta, prec);

  arb_neg(arb_mat_entry(M[2], 2, 0), sinalpha);
  arb_zero(arb_mat_entry(M[2], 2, 1));
  arb_set(arb_mat_entry(M[2], 2, 2), cosalpha);

  for (slong i = p->boundary; i < p->boundary + p->interior; i++)
  {
    /* We take s and t random with s in [0, 1) and s < 1 - t */
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

    for (slong vertex = 0; vertex < 3; vertex++)
    {
      /* Do a change of coordinates */
      arb_set(arb_mat_entry(xyz_matrix, 0, 0), xyz + 0);
      arb_set(arb_mat_entry(xyz_matrix, 1, 0), xyz + 1);
      arb_set(arb_mat_entry(xyz_matrix, 2, 0), xyz + 2);

      arb_mat_mul(xyz_matrix, M[vertex], xyz_matrix, prec);

      arb_set(xyz + 0, arb_mat_entry(xyz_matrix, 0, 0));
      arb_set(xyz + 1, arb_mat_entry(xyz_matrix, 1, 0));
      arb_set(xyz + 2, arb_mat_entry(xyz_matrix, 2, 0));

      cartesian_to_spherical(p->thetas[vertex] + i, p->phis[vertex] + i, xyz,
                             prec);

    }

  }

  arb_mat_clear(M[0]);
  arb_mat_clear(M[1]);
  arb_mat_clear(M[2]);

  _arb_vec_clear(xyz, 3);
  _arb_vec_clear(v2mv1, 3);
  _arb_vec_clear(v3mv1, 3);

  arb_clear(alpha);
  arb_clear(beta);
  arb_clear(cosalpha);
  arb_clear(cosbeta);
  arb_clear(sinalpha);
  arb_clear(sinbeta);
  arb_clear(s);
  arb_clear(t);
  arb_clear(tmp);
}
