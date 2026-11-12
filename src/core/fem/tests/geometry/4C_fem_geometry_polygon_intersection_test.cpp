// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fem_geometry_polygon_intersection.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace
{
  using namespace FourC;

  using Points = std::vector<std::array<double, 2>>;
  using Origin = Core::Geo::PolygonIntersectionPoint::Origin;

  const Points unit_square = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};

  // square with the given center and half edge length, rotated by angle
  Points square(double center_x, double center_y, double half_edge, double angle)
  {
    Points corners;
    for (int i = 0; i < 4; ++i)
    {
      const double phi = angle + (i + 0.5) * M_PI / 2.0;
      corners.push_back({center_x + std::sqrt(2.0) * half_edge * std::cos(phi),
          center_y + std::sqrt(2.0) * half_edge * std::sin(phi)});
    }
    return corners;
  }

  auto count(const std::vector<Core::Geo::PolygonIntersectionPoint>& intersection, Origin origin)
  {
    return std::ranges::count(intersection, origin, &Core::Geo::PolygonIntersectionPoint::origin);
  }

  TEST(IntersectConvexPolygons2D, RotatedSquaresIntersectInOctagon)
  {
    const Points rotated = square(0.5, 0.5, 0.5, M_PI / 4.0);
    const auto intersection = Core::Geo::intersect_convex_polygons_2d(unit_square, rotated, 1e-10);

    ASSERT_EQ(intersection.size(), 8);
    EXPECT_EQ(count(intersection, Origin::edge_intersection), 8);
    for (const auto& point : intersection)
    {
      // the intersection lies on the edges index_a of the unit square and index_b of the rotated
      // one
      const auto& start = rotated[point.index_b];
      const auto& end = rotated[(point.index_b + 1) % 4];
      const double distance_from_edge_b =
          ((end[0] - start[0]) * (point.coordinates[1] - start[1]) -
              (end[1] - start[1]) * (point.coordinates[0] - start[0])) /
          std::hypot(end[0] - start[0], end[1] - start[1]);
      EXPECT_NEAR(distance_from_edge_b, 0.0, 1e-14);
      const bool on_edge_a = (point.index_a == 0 and point.coordinates[1] == 0.0) or
                             (point.index_a == 1 and point.coordinates[0] == 1.0) or
                             (point.index_a == 2 and point.coordinates[1] == 1.0) or
                             (point.index_a == 3 and point.coordinates[0] == 0.0);
      EXPECT_TRUE(on_edge_a);
    }
  }

  TEST(IntersectConvexPolygons2D, PolygonInsideOtherPolygon)
  {
    const auto intersection =
        Core::Geo::intersect_convex_polygons_2d(unit_square, square(0.5, 0.5, 0.25, 0.3), 1e-10);

    ASSERT_EQ(intersection.size(), 4);
    EXPECT_EQ(count(intersection, Origin::corner_b), 4);
  }

  // Corners of a are kept before identical corners of b and intersections at the corners.
  TEST(IntersectConvexPolygons2D, IdenticalPolygonsGiveCornersOfA)
  {
    const auto intersection =
        Core::Geo::intersect_convex_polygons_2d(unit_square, unit_square, 1e-10);

    ASSERT_EQ(intersection.size(), 4);
    EXPECT_EQ(count(intersection, Origin::corner_a), 4);
  }

  // Almost identical polygons, e.g. from matching meshes, do not create intersections of almost
  // parallel edges somewhere in the middle of the edges.
  TEST(IntersectConvexPolygons2D, AlmostIdenticalPolygonsGiveCorners)
  {
    const auto intersection = Core::Geo::intersect_convex_polygons_2d(
        unit_square, square(0.5 + 1e-13, 0.5 - 1e-13, 0.5, 1e-13), 1e-10);

    ASSERT_EQ(intersection.size(), 4);
    EXPECT_EQ(count(intersection, Origin::corner_a), 4);
  }

  // A corner of b on an edge of a counts as inside and replaces the intersections at that corner.
  TEST(IntersectConvexPolygons2D, CornerOnEdge)
  {
    const Points triangle = {{0.5, 0.0}, {1.5, -1.0}, {1.5, 1.0}};
    const auto intersection = Core::Geo::intersect_convex_polygons_2d(unit_square, triangle, 1e-10);

    ASSERT_EQ(intersection.size(), 3);
    EXPECT_EQ(count(intersection, Origin::corner_b), 1);
    EXPECT_EQ(count(intersection, Origin::edge_intersection), 1);
    EXPECT_EQ(count(intersection, Origin::corner_a), 1);
  }

  TEST(IntersectConvexPolygons2D, NoAreaGivesEmptyIntersection)
  {
    const Points touching = {{1.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {1.0, 1.0}};
    const Points disjoint = {{3.0, 0.0}, {4.0, 0.0}, {4.0, 1.0}, {3.0, 1.0}};

    EXPECT_TRUE(Core::Geo::intersect_convex_polygons_2d(unit_square, touching, 1e-10).empty());
    EXPECT_TRUE(Core::Geo::intersect_convex_polygons_2d(unit_square, disjoint, 1e-10).empty());
  }

  TEST(IntersectConvexPolygons2D, NonConvexPolygonGivesEmptyIntersection)
  {
    const Points non_convex = {{0.0, 0.0}, {1.0, 0.0}, {0.5, 0.2}, {0.5, 1.0}};

    EXPECT_TRUE(Core::Geo::intersect_convex_polygons_2d(unit_square, non_convex, 1e-10).empty());
  }

  // Flat rectangle sharing its right edge with a tall rectangle (from a volume mortar test).
  TEST(IntersectConvexPolygons2D, SharedEdge)
  {
    const Points flat = {{0.0, 0.0}, {0.12264150943396501, 1.2325951644078309e-32},
        {0.12264150943396371, 0.024999999999999856}, {4.463523567848754e-15, 0.025000000000000022}};
    const Points tall = {{-0.36792452830188349, -0.024999999999999578},
        {0.1226415094339641, -0.025000000000000022}, {0.12264150943396455, 0.47499999999999998},
        {-0.36792452830188482, 0.47500000000000042}};
    const auto intersection = Core::Geo::intersect_convex_polygons_2d(flat, tall, 1e-12);

    ASSERT_EQ(intersection.size(), 4);
    EXPECT_EQ(count(intersection, Origin::corner_a), 4);
  }

  // The orientation of the polygons does not matter.
  TEST(IntersectConvexPolygons2D, ClockwisePolygons)
  {
    const Points clockwise_square = {{0.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}, {1.0, 0.0}};
    const auto intersection = Core::Geo::intersect_convex_polygons_2d(
        clockwise_square, square(0.5, 0.5, 0.5, M_PI / 4.0), 1e-10);

    EXPECT_EQ(intersection.size(), 8);
  }
}  // namespace
