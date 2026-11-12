// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_geometry_delaunay_triangulation.hpp"

#include <algorithm>
#include <cmath>

FOUR_C_NAMESPACE_OPEN

namespace
{
  // Check if a point lies inside the circumcircle of a triangle by more than the tolerance
  bool is_inside_circumcircle(const std::array<double, 2>& a, const std::array<double, 2>& b,
      const std::array<double, 2>& c, const std::array<double, 2>& point, double tolerance)
  {
    // circumcenter relative to a
    const double bx = b[0] - a[0];
    const double by = b[1] - a[1];
    const double cx = c[0] - a[0];
    const double cy = c[1] - a[1];
    const double d = 2.0 * (bx * cy - by * cx);
    const double center_x = (cy * (bx * bx + by * by) - by * (cx * cx + cy * cy)) / d;
    const double center_y = (bx * (cx * cx + cy * cy) - cx * (bx * bx + by * by)) / d;

    const double radius = std::hypot(center_x, center_y);
    const double distance = std::hypot(point[0] - a[0] - center_x, point[1] - a[1] - center_y);
    return distance < radius - tolerance;
  }
}  // namespace

std::vector<std::array<int, 3>> Core::Geo::delaunay_triangulation_convex_polygon_2d(
    std::span<const std::array<double, 2>> corners, double tolerance)
{
  std::vector<std::array<int, 3>> triangles;
  for (int i = 1; i + 1 < static_cast<int>(corners.size()); ++i) triangles.push_back({0, i, i + 1});

  // flip the shared edge of two triangles if the opposite corner of one lies inside the
  // circumcircle of the other, until no flip is necessary
  bool flipped = true;
  while (flipped)
  {
    flipped = false;
    for (std::size_t t = 0; t < triangles.size() and !flipped; ++t)
    {
      for (std::size_t u = t + 1; u < triangles.size() and !flipped; ++u)
      {
        const auto contains = [](const auto& range, int corner)
        { return std::ranges::find(range, corner) != std::ranges::end(range); };

        std::vector<int> shared;
        for (int corner : triangles[t])
          if (contains(triangles[u], corner)) shared.push_back(corner);
        if (shared.size() != 2) continue;

        const auto not_shared = [&](int corner) { return !contains(shared, corner); };
        const int opposite_t = *std::ranges::find_if(triangles[t], not_shared);
        const int opposite_u = *std::ranges::find_if(triangles[u], not_shared);

        if (is_inside_circumcircle(corners[triangles[t][0]], corners[triangles[t][1]],
                corners[triangles[t][2]], corners[opposite_u], tolerance))
        {
          triangles[t] = {opposite_t, shared[0], opposite_u};
          triangles[u] = {opposite_t, opposite_u, shared[1]};
          flipped = true;
        }
      }
    }
  }

  return triangles;
}

FOUR_C_NAMESPACE_CLOSE
