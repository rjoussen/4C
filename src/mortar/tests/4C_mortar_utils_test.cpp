// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_mortar_utils.hpp"

#include "4C_linalg_serialdensematrix.hpp"
#include "4C_mortar_coupling3d_classes.hpp"

#include <array>
#include <vector>

namespace
{
  using namespace FourC;

  // Points of a candidate clip polygon in the coordinates of the auxiliary plane
  using Points = std::vector<std::array<double, 2>>;

  // Sort the candidate points and return the resulting clip polygon
  Points sort_convex_hull(const Points& points, double tol = 1.0e-8)
  {
    Core::LinAlg::SerialDenseMatrix transformed(2, points.size());
    std::vector<Mortar::Vertex> candidates;
    for (std::size_t i = 0; i < points.size(); ++i)
    {
      transformed(0, i) = points[i][0];
      transformed(1, i) = points[i][1];
      candidates.emplace_back(std::vector<double>{points[i][0], points[i][1], 0.0},
          Mortar::Vertex::lineclip, std::vector<int>{}, nullptr, nullptr, false, false, nullptr,
          -1.0);
    }

    std::vector<Mortar::Vertex> polygon;
    Mortar::sort_convex_hull_points(false, transformed, candidates, polygon, tol);

    Points result;
    for (const auto& vertex : polygon) result.push_back({vertex.coord()[0], vertex.coord()[1]});
    return result;
  }

  double signed_area(const Points& polygon)
  {
    double area = 0.0;
    for (std::size_t i = 0; i < polygon.size(); ++i)
    {
      const auto& p = polygon[i];
      const auto& q = polygon[(i + 1) % polygon.size()];
      area += 0.5 * (p[0] * q[1] - q[0] * p[1]);
    }
    return area;
  }

  void expect_polygon(const Points& expected, const Points& actual)
  {
    ASSERT_EQ(expected.size(), actual.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
      EXPECT_NEAR(expected[i][0], actual[i][0], 1.0e-14) << "vertex " << i;
      EXPECT_NEAR(expected[i][1], actual[i][1], 1.0e-14) << "vertex " << i;
    }
  }

  // The clip polygon starts at the point with the smallest x-coordinate and runs clockwise.
  TEST(SortConvexHullPoints, ClockwiseFromSmallestX)
  {
    const Points square = {{1.0, 1.0}, {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
    const Points polygon = sort_convex_hull(square);

    expect_polygon({{0.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {1.0, 0.0}}, polygon);
    EXPECT_NEAR(signed_area(polygon), -1.0, 1.0e-14);
  }

  // A point almost on a polygon edge is removed, not the corner behind it.
  // Such points arise as intersections of slave and master edges that (almost) coincide.
  // Here, the point lies slightly inside the first edge, so its angle is slightly smaller and
  // it is sorted after the far corner of that edge.
  TEST(SortConvexHullPoints, NearlyCollinearPointOnFirstEdge)
  {
    const Points diamond_with_point_on_edge = {
        {0.0, 0.0}, {1.0, 1.0}, {0.5, 0.5 - 1.0e-10}, {2.0, 0.0}, {1.0, -1.0}};
    const Points polygon = sort_convex_hull(diamond_with_point_on_edge);

    expect_polygon({{0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}, {1.0, -1.0}}, polygon);
    EXPECT_NEAR(signed_area(polygon), -2.0, 1.0e-9);
  }

  // Same as above for the last edge, which closes the polygon.
  TEST(SortConvexHullPoints, NearlyCollinearPointOnLastEdge)
  {
    const Points diamond_with_point_on_edge = {
        {0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}, {1.0, -1.0}, {0.5, -0.5 + 1.0e-10}};
    const Points polygon = sort_convex_hull(diamond_with_point_on_edge);

    expect_polygon({{0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}, {1.0, -1.0}}, polygon);
    EXPECT_NEAR(signed_area(polygon), -2.0, 1.0e-9);
  }

  // Points clearly off an edge are kept as corners.
  TEST(SortConvexHullPoints, CornerBeyondToleranceIsKept)
  {
    const Points rectangle_with_corner = {
        {0.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {1.0, 0.0}, {0.5, 1.0 + 1.0e-6}};
    const Points polygon = sort_convex_hull(rectangle_with_corner);

    expect_polygon({{0.0, 0.0}, {0.0, 1.0}, {0.5, 1.0 + 1.0e-6}, {1.0, 1.0}, {1.0, 0.0}}, polygon);
  }

  // A real corner next to a short edge is kept. The point (0.92405832073, 4.779e-7) lies slightly
  // inside the short right edge, so it is not part of the convex hull and the corner
  // (0.92405832344, 0) must not be judged against it. This sliver occurs in
  // fsi_ves_mono_mtrss_ost_ost_slidrot.
  TEST(SortConvexHullPoints, CornerAtShortEdgeOfSliverIsKept)
  {
    const Points sliver = {{0.0, 0.0}, {9.24058323435245743e-01, 2.60208521396521064e-18},
        {9.24058326599087154e-01, 9.46784155580679410e-07},
        {2.47319341372934974e-08, 9.20895476335759461e-07},
        {9.24058320726020455e-01, 4.77886198589801137e-07}};
    const Points polygon = sort_convex_hull(sliver);

    expect_polygon({sliver[0], sliver[3], sliver[2], sliver[1]}, polygon);
  }

  // If several points form an almost straight line with their neighbors, the straightest one is
  // removed first. Here, both the point (8.77e-7, 0) and the corner (0, 0) of this tiny polygon
  // are below the tolerance, but only the former is (almost) on an edge. This polygon occurs in
  // elch_3D_tet4_s2i_butlervolmer_mortar_standard.
  TEST(SortConvexHullPoints, StraightestPointIsRemovedFirst)
  {
    const Points tiny = {{0.0, 0.0}, {8.765804024145833e-07, 2.6469779601696886e-23},
        {8.9339976370960088e-07, 2.2107403755386921e-06},
        {2.3087783999301522e-06, 2.3492923768170549e-08}};
    const Points polygon = sort_convex_hull(tiny, 2.0872355153175881e-12);

    expect_polygon({tiny[0], tiny[2], tiny[3]}, polygon);
  }
}  // namespace
