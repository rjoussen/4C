// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GEOMETRY_DELAUNAY_TRIANGULATION_HPP
#define FOUR_C_FEM_GEOMETRY_DELAUNAY_TRIANGULATION_HPP

#include "4C_config.hpp"

#include <array>
#include <span>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  /**
   * @brief Compute the Delaunay triangulation of a convex polygon in 2D.
   *
   * The triangulation starts as a fan from the first corner. Diagonals are flipped as long as a
   * corner lies inside the circumcircle of a neighboring triangle by more than the tolerance. If
   * several triangulations are Delaunay within the tolerance, e.g. for a rectangle, the one closest
   * to the fan is kept.
   *
   * @param corners corners of the polygon, consecutive along its boundary
   * @param tolerance distance by which a corner must lie inside a circumcircle to flip a diagonal
   * @return corner indices of the triangles
   */
  std::vector<std::array<int, 3>> delaunay_triangulation_convex_polygon_2d(
      std::span<const std::array<double, 2>> corners, double tolerance);
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
