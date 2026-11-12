// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_fem_geometry_delaunay_triangulation.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace
{
  using namespace FourC;

  using Points = std::vector<std::array<double, 2>>;
  using Triangles = std::vector<std::array<int, 3>>;

  bool has_edge(const Triangles& triangles, int a, int b)
  {
    return std::ranges::any_of(triangles,
        [&](const std::array<int, 3>& triangle)
        {
          return std::ranges::find(triangle, a) != triangle.end() and
                 std::ranges::find(triangle, b) != triangle.end();
        });
  }

  TEST(DelaunayTriangulationConvexPolygon2D, TriangleIsKept)
  {
    const Triangles triangles = Core::Geo::delaunay_triangulation_convex_polygon_2d(
        Points{{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}}, 1e-10);

    EXPECT_EQ(triangles, (Triangles{{0, 1, 2}}));
  }

  // Both diagonals of a square are Delaunay, so the fan from the first corner is kept.
  TEST(DelaunayTriangulationConvexPolygon2D, SquareKeepsFan)
  {
    const Triangles triangles = Core::Geo::delaunay_triangulation_convex_polygon_2d(
        Points{{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}, 1e-10);

    ASSERT_EQ(triangles.size(), 2);
    EXPECT_TRUE(has_edge(triangles, 0, 2));
  }

  // The long diagonal of a rhombus from the first corner is flipped to the short one.
  TEST(DelaunayTriangulationConvexPolygon2D, RhombusUsesShortDiagonal)
  {
    const Triangles triangles = Core::Geo::delaunay_triangulation_convex_polygon_2d(
        Points{{-2.0, 0.0}, {0.0, -1.0}, {2.0, 0.0}, {0.0, 1.0}}, 1e-10);

    ASSERT_EQ(triangles.size(), 2);
    EXPECT_TRUE(has_edge(triangles, 1, 3));
    EXPECT_FALSE(has_edge(triangles, 0, 2));
  }

  // No corner of an irregular octagon lies inside the circumcircle of a triangle.
  TEST(DelaunayTriangulationConvexPolygon2D, OctagonIsDelaunay)
  {
    Points corners;
    const std::array<double, 8> angles = {0.0, 0.3, 1.2, 1.5, 2.9, 3.6, 4.4, 5.9};
    const std::array<double, 8> radii = {1.0, 1.1, 0.9, 1.0, 1.2, 0.8, 1.0, 1.1};
    for (int i = 0; i < 8; ++i)
      corners.push_back({radii[i] * std::cos(angles[i]), radii[i] * std::sin(angles[i])});

    const Triangles triangles = Core::Geo::delaunay_triangulation_convex_polygon_2d(corners, 0.0);

    ASSERT_EQ(triangles.size(), 6);
    for (const auto& [a, b, c] : triangles)
    {
      const double bx = corners[b][0] - corners[a][0];
      const double by = corners[b][1] - corners[a][1];
      const double cx = corners[c][0] - corners[a][0];
      const double cy = corners[c][1] - corners[a][1];
      const double d = 2.0 * (bx * cy - by * cx);
      const double center_x = (cy * (bx * bx + by * by) - by * (cx * cx + cy * cy)) / d;
      const double center_y = (bx * (cx * cx + cy * cy) - cx * (bx * bx + by * by)) / d;
      const double radius = std::hypot(center_x, center_y);

      for (const auto& corner : corners)
      {
        const double distance =
            std::hypot(corner[0] - corners[a][0] - center_x, corner[1] - corners[a][1] - center_y);
        EXPECT_GE(distance, radius - 1e-12);
      }
    }
  }
}  // namespace
