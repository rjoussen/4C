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
  Points sort_convex_hull(const Points& points)
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
    double tol = 1.0e-8;
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

  // A point within the tolerance of a polygon edge is removed, not the corner behind it.
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

  // Points closer to an edge than the tolerance are removed, independent of the edge length. This
  // is consistent with the inside/outside checks of the polygon clipping, which also use the
  // tolerance as a distance.
  TEST(SortConvexHullPoints, CollinearToleranceIsADistance)
  {
    const Points long_rectangle_with_point_near_edge = {
        {0.0, 0.0}, {0.0, 1.0}, {10.0, 1.0}, {10.0, 0.0}, {5.0, 1.0 + 5.0e-9}};
    const Points polygon = sort_convex_hull(long_rectangle_with_point_near_edge);

    expect_polygon({{0.0, 0.0}, {0.0, 1.0}, {10.0, 1.0}, {10.0, 0.0}}, polygon);
  }

  // Points farther from an edge than the tolerance are kept as corners.
  TEST(SortConvexHullPoints, CornerBeyondToleranceIsKept)
  {
    const Points rectangle_with_corner = {
        {0.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {1.0, 0.0}, {0.5, 1.0 + 1.0e-6}};
    const Points polygon = sort_convex_hull(rectangle_with_corner);

    expect_polygon({{0.0, 0.0}, {0.0, 1.0}, {0.5, 1.0 + 1.0e-6}, {1.0, 1.0}, {1.0, 0.0}}, polygon);
  }
}  // namespace
