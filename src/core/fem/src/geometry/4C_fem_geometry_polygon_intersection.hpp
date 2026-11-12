// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GEOMETRY_POLYGON_INTERSECTION_HPP
#define FOUR_C_FEM_GEOMETRY_POLYGON_INTERSECTION_HPP

#include "4C_config.hpp"

#include <array>
#include <span>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  /**
   * @brief A corner of the intersection of two polygons, together with its origin.
   */
  struct PolygonIntersectionPoint
  {
    enum class Origin
    {
      corner_a,          ///< corner index_a of polygon a
      corner_b,          ///< corner index_b of polygon b
      edge_intersection  ///< intersection of edge index_a of a and edge index_b of b
    };

    Origin origin;

    /// corner or edge of polygon a (edge i connects corners i and i+1), -1 if not involved
    int index_a;

    /// corner or edge of polygon b (edge i connects corners i and i+1), -1 if not involved
    int index_b;

    std::array<double, 2> coordinates;
  };

  /**
   * @brief Compute the intersection of two convex polygons in 2D.
   *
   * The corners of the intersection are corners of one polygon inside the other polygon and
   * intersections of edges. Corners closer than the tolerance to the other polygon count as
   * inside. Points closer than the tolerance to each other are merged, keeping corners of a before
   * corners of b before edge intersections. Edges that deviate less than the tolerance from being
   * parallel along their length do not create intersections, since their corners already describe
   * the intersection.
   *
   * @param polygon_a corners of the first polygon in either orientation
   * @param polygon_b corners of the second polygon in either orientation
   * @param tolerance distance below which points are treated as identical
   * @return corners of the intersection, consecutive along its boundary (empty if one of the
   * polygons is not convex or if the intersection has no area)
   */
  std::vector<PolygonIntersectionPoint> intersect_convex_polygons_2d(
      std::span<const std::array<double, 2>> polygon_a,
      std::span<const std::array<double, 2>> polygon_b, double tolerance);
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
