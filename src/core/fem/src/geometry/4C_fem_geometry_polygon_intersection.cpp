// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_geometry_polygon_intersection.hpp"

#include "4C_fem_geometry_convex_hull.hpp"

#include <boost/geometry/algorithms/correct.hpp>
#include <boost/geometry/algorithms/distance.hpp>
#include <boost/geometry/algorithms/intersection.hpp>
#include <boost/geometry/algorithms/is_convex.hpp>
#include <boost/geometry/algorithms/within.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/segment.hpp>

#include <algorithm>
#include <cmath>

FOUR_C_NAMESPACE_OPEN

namespace
{
  using Point = boost::geometry::model::d2::point_xy<double>;
  using Polygon = boost::geometry::model::polygon<Point>;
  using Linestring = boost::geometry::model::linestring<Point>;
  using Segment = boost::geometry::model::segment<Point>;

  Point to_point(const std::array<double, 2>& coordinates)
  {
    return {coordinates[0], coordinates[1]};
  }

  Polygon to_polygon(std::span<const std::array<double, 2>> corners)
  {
    Polygon polygon;
    for (const auto& corner : corners) polygon.outer().push_back(to_point(corner));
    boost::geometry::correct(polygon);
    return polygon;
  }

  // Check if a point is inside a polygon or closer than the tolerance to its boundary. The points
  // on the boundary are found by distance instead of boost::geometry::covered_by, which counts
  // points on the extension of almost vertical edges as boundary points.
  bool is_inside(const Point& point, const Polygon& polygon, double tolerance)
  {
    const Linestring boundary(polygon.outer().begin(), polygon.outer().end());
    return boost::geometry::within(point, polygon) or
           boost::geometry::distance(point, boundary) <= tolerance;
  }
}  // namespace

std::vector<Core::Geo::PolygonIntersectionPoint> Core::Geo::intersect_convex_polygons_2d(
    std::span<const std::array<double, 2>> polygon_a,
    std::span<const std::array<double, 2>> polygon_b, double tolerance)
{
  using Origin = PolygonIntersectionPoint::Origin;

  const Polygon boost_polygon_a = to_polygon(polygon_a);
  const Polygon boost_polygon_b = to_polygon(polygon_b);
  if (!boost::geometry::is_convex(boost_polygon_a) or !boost::geometry::is_convex(boost_polygon_b))
    return {};

  // candidates closer than the tolerance to a previous candidate are merged into it
  std::vector<PolygonIntersectionPoint> candidates;
  const auto add_candidate = [&](const PolygonIntersectionPoint& candidate)
  {
    for (const PolygonIntersectionPoint& previous : candidates)
    {
      if (std::hypot(candidate.coordinates[0] - previous.coordinates[0],
              candidate.coordinates[1] - previous.coordinates[1]) <= tolerance)
        return;
    }
    candidates.push_back(candidate);
  };

  for (int i = 0; i < static_cast<int>(polygon_a.size()); ++i)
  {
    if (is_inside(to_point(polygon_a[i]), boost_polygon_b, tolerance))
      add_candidate(
          {.origin = Origin::corner_a, .index_a = i, .index_b = -1, .coordinates = polygon_a[i]});
  }

  for (int j = 0; j < static_cast<int>(polygon_b.size()); ++j)
  {
    if (is_inside(to_point(polygon_b[j]), boost_polygon_a, tolerance))
      add_candidate(
          {.origin = Origin::corner_b, .index_a = -1, .index_b = j, .coordinates = polygon_b[j]});
  }

  for (int i = 0; i < static_cast<int>(polygon_a.size()); ++i)
  {
    const auto& start_a = polygon_a[i];
    const auto& end_a = polygon_a[(i + 1) % polygon_a.size()];

    for (int j = 0; j < static_cast<int>(polygon_b.size()); ++j)
    {
      const auto& start_b = polygon_b[j];
      const auto& end_b = polygon_b[(j + 1) % polygon_b.size()];

      // the distance between almost parallel edges changes by less than the tolerance along the
      // shorter edge
      const double cross_product = (end_a[0] - start_a[0]) * (end_b[1] - start_b[1]) -
                                   (end_a[1] - start_a[1]) * (end_b[0] - start_b[0]);
      const double longer_edge = std::max(std::hypot(end_a[0] - start_a[0], end_a[1] - start_a[1]),
          std::hypot(end_b[0] - start_b[0], end_b[1] - start_b[1]));
      if (std::abs(cross_product) <= tolerance * longer_edge) continue;

      std::vector<Point> intersection;
      boost::geometry::intersection(Segment(to_point(start_a), to_point(end_a)),
          Segment(to_point(start_b), to_point(end_b)), intersection);
      if (intersection.size() == 1)
        add_candidate({.origin = Origin::edge_intersection,
            .index_a = i,
            .index_b = j,
            .coordinates = {intersection[0].x(), intersection[0].y()}});
    }
  }

  std::vector<std::array<double, 2>> candidate_coordinates;
  for (const PolygonIntersectionPoint& candidate : candidates)
    candidate_coordinates.push_back(candidate.coordinates);

  std::vector<PolygonIntersectionPoint> intersection;
  for (int index : convex_hull_2d(candidate_coordinates)) intersection.push_back(candidates[index]);
  return intersection;
}

FOUR_C_NAMESPACE_CLOSE
